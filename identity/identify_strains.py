#!/usr/bin/env python3
"""Report closest reference genomes for local assemblies by ANI.

The diamond/MEGAN call in each strain folder is useful search-scope evidence but is
not a taxonomic verdict. This script returns ANI hits and aligned fractions so the
result can be interpreted manually; a closest reference below the species threshold
must not be written back as a species identification.

Per genus the reference set is one RefSeq genome per species:

1. NCBI's ``reference_only`` set, which is curated at roughly one assembly per
   species. Sparse for under-sequenced genera (Gluconobacter returns one).
2. When that yields fewer than ``--min-species`` species, a paginated scan of the
   genus grouped by species binomial, keeping the best assembly per species
   (reference/representative > complete genome > largest contig N50).

Taxon *names* are not used for the assembly lookup: NCBI resolves "Bacillus" to
three assemblies because the genus has been split, whereas tax id 1386 resolves
to 133. Genus tax ids are therefore resolved at runtime and validated by rank and
lineage before use.

The search is **exhaustive by default**: every species in the genus, and every
species in the family on escalation. References are sketched, indexed and queried
in batches so that peak memory is set by ``--batch-size`` rather than by how
well-sequenced the clade happens to be. The earlier 80-reference cap and size/GC
window were memory guards, not taxonomy, and both could hide the right answer --
the cap left 123 of 203 Enterobacteriaceae species unsearched for ``P2_G9``. They
remain available behind ``--prefilter`` and ``--max-references`` for a quick pass.

The ``call`` column is a reading aid, not an automated gate. >=95% ANI is labelled a
species match, >=99.5% with >=80% alignment is labelled a close strain-level match,
and 90-95% is only a genus-range hint. The hit table remains the evidence.

Writes results/identity/strain_identity.tsv and results/identity/ani_all_hits.tsv.
"""

from __future__ import annotations

import argparse
import csv
import http.client
import io
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from collections.abc import Iterator
from pathlib import Path

import pyfastani

# prepare_candidates.py copies the assembly, the QUAST report and the MEGAN taxonomy for
# every panel member, so once a strain is in the panel this script no longer needs the
# drive. The staged copies are the default; use --genome-root explicitly for another
# cohort.
ROOT = Path(__file__).resolve().parents[1]
LOCAL_GENOME_ROOT = ROOT / "candidate_genomes"

BACTERIA_TAX_ID = 2

SPECIES_ANI = 95.0
STRAIN_ANI = 99.5
GENUS_ANI = 90.0
MIN_ALIGNED_FRACTION = 0.35
STRAIN_MIN_ALIGNED_FRACTION = 0.80

USER_AGENT = "glycerol-lactate-ani/1.0"
REPORT_URL = (
    "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{tax_id}/dataset_report"
    "?page_size={page_size}{extra}"
)
TAXONOMY_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/{query}"
DOWNLOAD_URL = (
    "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}"
    "/download?include_annotation_type=GENOME_FASTA"
)
FASTA_SUFFIXES = {".fa", ".fasta", ".fna"}


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


# NCBI asks for no more than 3 requests/second without an API key. Reference
# downloads are one request per genome, so a genus scan is hundreds of calls —
# throttle every request through one gate and back off exponentially on failure.
MIN_REQUEST_INTERVAL_SECONDS = 0.4
_last_request_time = 0.0


def throttle() -> None:
    global _last_request_time
    elapsed = time.monotonic() - _last_request_time
    if elapsed < MIN_REQUEST_INTERVAL_SECONDS:
        time.sleep(MIN_REQUEST_INTERVAL_SECONDS - elapsed)
    _last_request_time = time.monotonic()


def request_bytes(url: str, attempts: int = 4) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    for attempt in range(1, attempts + 1):
        throttle()
        try:
            with urllib.request.urlopen(request, timeout=180) as response:
                return response.read()
        except urllib.error.HTTPError as error:
            # 429/5xx are worth retrying with a longer pause; 4xx are not.
            if error.code not in (429, 500, 502, 503, 504) or attempt == attempts:
                raise RuntimeError(f"Could not download {url}: {error}") from error
            log(f"    HTTP {error.code}, backing off {2 ** attempt}s")
            time.sleep(2**attempt)
        except (urllib.error.URLError, TimeoutError, http.client.HTTPException,
                ConnectionError, OSError) as error:
            # IncompleteRead in particular: a chunked response truncated mid-transfer.
            # It is not an HTTPError, so it slipped past the retry above and killed a
            # 436-genome family search on its 300th download.
            if attempt == attempts:
                raise RuntimeError(f"Could not download {url}: {error}") from error
            time.sleep(2**attempt)
    raise AssertionError("unreachable")


def taxonomy_nodes(query: str) -> list[dict]:
    payload = json.loads(request_bytes(TAXONOMY_URL.format(query=urllib.parse.quote(str(query)))))
    return [node.get("taxonomy", {}) for node in payload.get("taxonomy_nodes", []) if node.get("taxonomy")]


def resolve_genus_tax_id(genus: str) -> int | None:
    """Genus name to tax id, validated by rank and bacterial lineage.

    Thin wrapper over resolve_taxon_tax_id, kept because the genus-only case is the
    default and the reasoning below is what makes it safe.

    Hand-pinning tax ids is error-prone in both directions: 442 looks like the
    Gluconobacter genus but is *G. oxydans*, so a genus reference set silently
    collapsed to one species. Names alone are also unsafe — "Bacillus" resolves to
    both the bacterial genus (1386) and a stick-insect genus (55087), and querying
    the dataset endpoint by name returned three assemblies instead of 133. So
    resolve by name, then require rank GENUS and Bacteria in the lineage.
    """
    return resolve_taxon_tax_id(genus, ("GENUS",))


def resolve_taxon_tax_id(name: str, ranks: tuple[str, ...]) -> int | None:
    """Resolve a taxon name at any of ``ranks``, with the same safety checks as above.

    Needed because the MEGAN call for some assemblies lands on a *strain* identifier
    that the taxonomy database carries at genus level -- ``P2_G9`` resolves to
    "BIGb0383", which is not a genus and has no reference genomes. Its deeper lineage is
    Enterobacteriaceae at 99.95%, so the useful move is to name the family explicitly and
    run ANI against that. Kept separate from resolve_genus_tax_id so the strict
    genus-only path stays the default and overrides have to be asked for.
    """
    for node in taxonomy_nodes(name):
        if (
            node.get("rank") in ranks
            and str(node.get("organism_name", "")).casefold() == name.casefold()
            and BACTERIA_TAX_ID in (node.get("lineage") or [])
        ):
            return int(node["tax_id"])
    return None


def resolve_family_tax_id(genus_tax_id: int) -> tuple[int, str] | None:
    """The family containing a genus, for escalation when no species matches."""
    nodes = taxonomy_nodes(str(genus_tax_id))
    if not nodes:
        return None
    lineage = [t for t in (nodes[0].get("lineage") or []) if t not in (1, 131567, BACTERIA_TAX_ID)]
    if not lineage:
        return None
    for ancestor in taxonomy_nodes(",".join(str(t) for t in reversed(lineage))):
        if ancestor.get("rank") == "FAMILY":
            return int(ancestor["tax_id"]), str(ancestor.get("organism_name", ""))
    return None


def merge_rows(
    path: Path,
    new_rows: list[dict],
    key,
    replace_scope=None,
    sort_key=None,
) -> None:
    """Write ``new_rows`` into ``path``, keeping rows for genomes not in this run.

    The script has to be invoked one genome per process -- that is what keeps peak memory
    to a single genus reference set instead of all of them at once. Writing in "w" mode
    made that pattern quietly destructive: each run replaced the whole table, so
    identifying one new strain erased the other eighteen. Rows are keyed and merged
    instead, and any genome touched by this run has its old rows dropped first so a
    re-identification does not leave stale hits behind.
    """
    if not new_rows:
        log(f"nothing to merge into {path.name}; leaving it untouched")
        return
    existing: list[dict] = []
    if path.exists():
        with path.open(newline="", encoding="utf-8") as handle:
            existing = list(csv.DictReader(handle, delimiter="\t"))
    scope = {(replace_scope or key)(r) for r in new_rows}
    kept = [r for r in existing if (replace_scope or key)(r) not in scope]
    merged = {key(r): r for r in kept}
    for row in new_rows:
        merged[key(row)] = row
    rows = list(merged.values())
    rows.sort(key=sort_key or key)
    fieldnames = list(new_rows[0].keys())
    for row in rows:                       # tolerate a schema that grew between runs
        for extra in row:
            if extra not in fieldnames:
                fieldnames.append(extra)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", restval="")
        writer.writeheader()
        writer.writerows(rows)


def fasta_contigs(path: Path) -> Iterator[bytes]:
    """Yield uppercase contigs without loading the whole file."""
    sequence = bytearray()
    with path.open("rb") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(b">"):
                if sequence:
                    yield bytes(sequence).upper()
                    sequence.clear()
            else:
                sequence.extend(line)
    if sequence:
        yield bytes(sequence).upper()


def assembly_priority(relative: Path) -> tuple[int, int, str]:
    """Rank duplicate assemblies within a strain folder; lower is preferred."""
    parts = tuple(part.lower() for part in relative.parts)
    filename = relative.name.lower()
    if len(relative.parts) == 2:
        rank = 0
    elif "spades_assembled" in parts and "pilon" in filename:
        rank = 1
    elif "spades_assembled" in parts:
        rank = 2
    elif "prokka" in parts:
        rank = 3
    else:
        rank = 4
    return rank, 0 if relative.suffix.lower() == ".fna" else 1, relative.as_posix()


def select_assembly(strain_dir: Path) -> Path | None:
    candidates = [
        path
        for path in strain_dir.rglob("*")
        if path.is_file()
        and path.suffix.lower() in FASTA_SUFFIXES
        and not path.name.startswith("._")
        and not any(part.startswith(".") for part in path.relative_to(strain_dir).parts)
    ]
    if not candidates:
        return None
    return min(candidates, key=lambda p: assembly_priority(Path(strain_dir.name) / p.relative_to(strain_dir)))


def megan_lineage(strain_dir: Path) -> list[tuple[str, float, int]]:
    """Best-supported taxon per rank row: (name, percent, number of competing names).

    The percentages are normalised **within each rank row**, so a row holding a single
    name always reads 100% no matter how few reads support it. That is why the third
    return value matters: a 100% call from a one-entry row is not evidence, and treating
    it as such is exactly how `P1_B3` came to be labelled "Aneurinibacillus, 100% genus
    support" while its own family row said Bacillaceae at 99.98%.
    """
    for base in (strain_dir / "diamond_megan", strain_dir):
        taxa = base / "unique_taxa.txt"
        pct = base / "unique_taxa_percentage.txt"
        if not taxa.exists():
            continue
        tax_lines = [x.strip().strip("|") for x in taxa.read_text(errors="replace").splitlines() if x.strip()]
        pct_lines = ([x.strip().strip("|") for x in pct.read_text(errors="replace").splitlines() if x.strip()]
                     if pct.exists() else [])
        lineage: list[tuple[str, float, int]] = []
        for index, line in enumerate(tax_lines):
            names = [n for n in line.split("|") if n]
            if not names:
                continue
            percentages: list[float] = []
            if index < len(pct_lines):
                try:
                    percentages = [float(x) for x in pct_lines[index].split("|") if x]
                except ValueError:
                    percentages = []
            if percentages and len(percentages) == len(names):
                name, percent = max(zip(names, percentages), key=lambda pair: pair[1])
            else:
                name, percent = names[0], 0.0
            lineage.append((name, percent, len(names)))
        if lineage:
            return lineage
    return []


def megan_genus(strain_dir: Path) -> tuple[str, str]:
    """Genus and its support percentage from the diamond/MEGAN output."""
    lineage = megan_lineage(strain_dir)
    if not lineage:
        return "", ""
    name, percent, _ = lineage[-1]
    return name, f"{percent:.2f}" if percent else ""


def megan_family(strain_dir: Path) -> tuple[str, float] | None:
    """The family MEGAN itself supports, which is not always the genus's family.

    Escalating from the genus assumes the genus call is right at family level. When it is
    not, escalation searches the wrong family and returns nothing -- `P1_B3` escalated from
    *Aneurinibacillus* into all 436 Paenibacillaceae species and matched **none** of them,
    while ANI places it in Bacillaceae at 99.47% (*Priestia aryabhattai*). Its own MEGAN
    output said so all along: family Bacillaceae 99.98% against Aneurinibacillaceae 0.02%.

    Read the family row directly and let it compete with the genus-derived family.
    """
    lineage = megan_lineage(strain_dir)
    for name, percent, competitors in lineage:
        if name.endswith("aceae") and competitors > 1 and percent > 0:
            return name, percent
    return None


def quast_metrics(strain_dir: Path) -> dict[str, str]:
    """Assembly QC from the existing QUAST output — an independent identity check.

    Total length and GC% are cheap corroboration: an assembly called Lactococcus
    that is 6 Mbp at 62% GC is not a Lactococcus regardless of what ANI says.
    """
    report = strain_dir / "quast_output" / "report.tsv"
    wanted = {
        "Total length": "assembly_bp",
        "GC (%)": "assembly_gc_percent",
        "# contigs": "contigs",
        "N50": "n50",
    }
    out = {value: "" for value in wanted.values()}
    if not report.exists():
        return out
    for line in report.read_text(errors="replace").splitlines():
        parts = line.split("\t")
        if len(parts) >= 2 and parts[0] in wanted:
            out[wanted[parts[0]]] = parts[1]
    return out


def fetch_reports(tax_id: int, extra: str, page_size: int = 1000, max_reports: int = 4000) -> list[dict]:
    reports: list[dict] = []
    token = ""
    while len(reports) < max_reports:
        url = REPORT_URL.format(tax_id=tax_id, page_size=page_size, extra=extra)
        if token:
            url += f"&page_token={token}"
        payload = json.loads(request_bytes(url))
        batch = payload.get("reports", [])
        if not batch:
            break
        reports.extend(batch)
        token = payload.get("next_page_token", "")
        if not token:
            break
    return reports


def species_of(report: dict) -> str:
    name = str(report.get("organism", {}).get("organism_name", "")).split()
    return " ".join(name[:2]) if len(name) >= 2 else " ".join(name)


def report_rank(report: dict) -> tuple[int, int, int]:
    info = report.get("assembly_info", {})
    category = str(info.get("refseq_category", "")).lower()
    category_rank = 0 if "reference" in category else (1 if "representative" in category else 2)
    level_rank = 0 if info.get("assembly_level") == "Complete Genome" else 1
    n50 = int(report.get("assembly_stats", {}).get("contig_n50", 0) or 0)
    return category_rank, level_rank, -n50


def build_reference_set(tax_id: int, genus: str, min_species: int) -> list[dict]:
    """One RefSeq assembly per species: curated set first, species-grouped scan if sparse."""
    curated = fetch_reports(tax_id, "&filters.reference_only=true")
    by_species: dict[str, dict] = {}
    for report in curated:
        name = species_of(report)
        if name and (name not in by_species or report_rank(report) < report_rank(by_species[name])):
            by_species[name] = report
    log(f"  {genus}: curated reference set covers {len(by_species)} species")

    if len(by_species) < min_species:
        log(f"  {genus}: sparse, scanning genus grouped by species")
        scanned = fetch_reports(tax_id, "&filters.assembly_source=refseq")
        for report in scanned:
            name = species_of(report)
            if not name:
                continue
            if name not in by_species or report_rank(report) < report_rank(by_species[name]):
                by_species[name] = report
        log(f"  {genus}: expanded to {len(by_species)} species")

    references = []
    for name, report in sorted(by_species.items()):
        accession = str(report.get("accession", ""))
        paired = str(report.get("paired_accession", ""))
        # Prefer the RefSeq member of a GenBank/RefSeq pair.
        if accession.startswith("GCA_") and paired.startswith("GCF_"):
            accession = paired
        stats = report.get("assembly_stats", {}) or {}
        references.append({
            "accession": accession,
            "species": name,
            "organism": str(report.get("organism", {}).get("organism_name", "")),
            "level": str(report.get("assembly_info", {}).get("assembly_level", "")),
            "category": str(report.get("assembly_info", {}).get("refseq_category", "")),
            "length": int(stats.get("total_sequence_length", 0) or 0),
            "gc_percent": float(stats.get("gc_percent", 0) or 0),
        })
    return references


def prefilter_by_size_and_gc(
    references: list[dict],
    query_length: int,
    query_gc: float,
    size_tolerance: float = 0.20,
    gc_tolerance: float = 4.0,
) -> list[dict]:
    """Drop references whose genome size or GC rules out a species-level match.

    Sketching every species in a large genus is the dominant cost of this script:
    indexing the 329-species Paenibacillus set means hashing 2.0 GB and holding ~1.4 GB
    resident, which repeatedly dominated the run and starved the machine — all to place
    a single query.

    Two organisms at >=95% ANI have similar genome size and GC, so references far
    outside the query's range cannot be the answer and need not be sketched. This is a
    conservative prefilter: the tolerances are wide relative to real intra-species
    variation, and if it removes everything the full set is used instead so a
    surprising genome is never silently excluded.
    """
    if not query_length or not query_gc:
        return references
    kept = [
        r for r in references
        if r.get("length") and r.get("gc_percent")
        and abs(r["length"] - query_length) <= size_tolerance * query_length
        and abs(r["gc_percent"] - query_gc) <= gc_tolerance
    ]
    if not kept:
        log("    prefilter removed every reference; falling back to the full set")
        return references
    return kept


def cap_references(references: list[dict], query_length: int, query_gc: float, limit: int) -> list[dict]:
    """Keep at most ``limit`` references, nearest in size and GC to the query.

    A hard cap makes peak memory predictable. Escalating a sparse genus to its family
    is otherwise unbounded — Aneurinibacillus has 9 species but Paenibacillaceae has
    436, so a failed genus match quietly turns into a multi-gigabyte download and
    index on a machine that has already been asked for 2 GB.
    """
    if len(references) <= limit or not query_length or not query_gc:
        return references[:limit]

    def distance(reference: dict) -> float:
        size_term = abs(reference.get("length", 0) - query_length) / max(query_length, 1)
        gc_term = abs(reference.get("gc_percent", 0.0) - query_gc) / 10.0
        return size_term + gc_term

    return sorted(references, key=distance)[:limit]


def download_reference(accession: str, reference_dir: Path) -> Path | None:
    destination = reference_dir / f"{re.sub(r'[^A-Za-z0-9._-]+', '_', accession)}.fna"
    if destination.exists() and destination.stat().st_size > 0:
        return destination
    try:
        archive = request_bytes(DOWNLOAD_URL.format(accession=accession))
        with zipfile.ZipFile(io.BytesIO(archive)) as package:
            members = [n for n in package.namelist() if n.endswith("_genomic.fna") and "/data/" in n]
            if not members:
                return None
            with package.open(members[0]) as source, destination.open("wb") as target:
                while chunk := source.read(1024 * 1024):
                    target.write(chunk)
    except (RuntimeError, zipfile.BadZipFile, OSError) as error:
        log(f"    WARNING: {accession} unavailable: {error}")
        destination.unlink(missing_ok=True)      # never leave a truncated file to be cached
        return None
    return destination


def classify(ani: float, aligned_fraction: float) -> str:
    if aligned_fraction < MIN_ALIGNED_FRACTION:
        return "low_alignment_unreliable"
    if ani >= STRAIN_ANI and aligned_fraction >= STRAIN_MIN_ALIGNED_FRACTION:
        return "strain_level_match"
    if ani >= SPECIES_ANI:
        return "species_match"
    if ani >= GENUS_ANI:
        return "genus_match_species_unresolved"
    return "no_species_match"


def search_references(
    label: str,
    tax_id: int,
    assembly: Path,
    options: argparse.Namespace,
    query_length: int = 0,
    query_gc: float = 0.0,
) -> list[tuple[float, float, str, dict]]:
    """ANI of ``assembly`` against every species under ``tax_id``, best first.

    References are sketched, indexed and queried in **batches**, and each batch's index is
    discarded before the next is built. That is what makes an exhaustive search affordable:
    the earlier version put the whole reference set into a single index, so memory scaled
    with the size of the genus or family and a 436-species family exhausted a 16 GB
    machine. Two guards were bolted on to survive that -- a size/GC prefilter and an
    80-reference cap -- and both are answers to a memory problem rather than to a
    taxonomic one.

    Both are now off by default, because both could silently hide the right answer. The
    cap dropped 123 of 203 Enterobacteriaceae species for `P2_G9`, so "no species match"
    rested on a search that never looked at most of the family. The GC/size window is
    worse in principle: it assumes the correct species resembles the query in exactly the
    properties that vary within a genus, which is precisely the assumption a
    reclassified or novel organism breaks.

    Peak memory is now set by ``--batch-size`` rather than by how well-sequenced the
    clade happens to be.
    """
    reference_dir = options.reference_dir / re.sub(r"[^A-Za-z0-9._-]+", "_", label)
    reference_dir.mkdir(parents=True, exist_ok=True)
    references = build_reference_set(tax_id, label, options.min_species)
    if options.prefilter:
        before = len(references)
        references = prefilter_by_size_and_gc(references, query_length, query_gc)
        if options.max_references:
            references = cap_references(references, query_length, query_gc, options.max_references)
        log(f"  {label}: size/GC prefilter kept {len(references)}/{before} references "
            f"(opt-in; off by default)")
    else:
        log(f"  {label}: searching all {len(references)} species, unfiltered and uncapped")

    hits: list[tuple[float, float, str, dict]] = []
    searched = 0
    skipped: list[str] = []
    for start in range(0, len(references), options.batch_size):
        batch = references[start:start + options.batch_size]
        sketch = pyfastani.Sketch()
        metadata: dict[str, dict] = {}
        for reference in batch:
            path = download_reference(reference["accession"], reference_dir)
            if path is None:
                skipped.append(reference["accession"])
                continue
            sketch.add_draft(reference["accession"], fasta_contigs(path))
            metadata[reference["accession"]] = reference
        if not metadata:
            continue
        mapper = sketch.index()
        hits.extend(best_hits(mapper, metadata, assembly, options.threads))
        searched += len(metadata)
        del mapper, sketch, metadata
        log(f"  {label}: searched {searched}/{len(references)}")
    if not searched:
        log(f"  {label}: no usable references")
    if skipped:
        # Say what was missed. A search that quietly covers less than it claims is how
        # "no species match" starts meaning "we did not look" -- which is the whole
        # reason the reference cap was removed.
        log(f"  {label}: WARNING {len(skipped)} reference(s) could not be downloaded and "
            f"were NOT searched: {', '.join(skipped[:5])}"
            f"{' ...' if len(skipped) > 5 else ''}")
    return sorted(hits, key=lambda item: item[0], reverse=True)


def best_hits(
    mapper: pyfastani.Mapper,
    metadata: dict[str, dict],
    assembly: Path,
    threads: int,
) -> list[tuple[float, float, str, dict]]:
    """(ani, aligned_fraction, reference_accession, reference_metadata), best first."""
    results = []
    for hit in mapper.query_draft(fasta_contigs(assembly), threads=threads):
        fraction = hit.matches / hit.fragments if hit.fragments else 0.0
        results.append((float(hit.identity), fraction, hit.name, metadata.get(hit.name, {})))
    return sorted(results, key=lambda item: item[0], reverse=True)


def parse_args() -> argparse.Namespace:
    root = ROOT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome-root", type=Path, default=LOCAL_GENOME_ROOT)
    parser.add_argument("--genome", action="append", help="Strain folder to identify (repeatable).")
    parser.add_argument("--manifest", type=Path,
                        default=root / "candidate_genomes" / "manifest_pool_v5.tsv",
                        help="Default genome list when --genome is not given.")
    parser.add_argument("--reference-dir", type=Path, default=root / "references")
    parser.add_argument("--output-dir", type=Path, default=root / "results" / "identity")
    parser.add_argument("--fallback-family", default=None,
                        help="Family reference set to search for genomes that carry no "
                             "MEGAN taxonomy at all.")
    parser.add_argument("--min-species", type=int, default=5,
                        help="Below this many species in the curated set, scan the genus.")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--batch-size", type=int, default=40,
                        help="References sketched and indexed at a time. Sets peak memory; "
                             "the search itself is exhaustive regardless.")
    parser.add_argument("--max-references", type=int, default=0,
                        help="Optional cap on references per genus/family, applied only with "
                             "--prefilter. 0 (default) means no cap: the search is exhaustive.")
    parser.add_argument("--prefilter", action="store_true",
                        help="Opt in to the old size/GC window. Off by default because it "
                             "assumes the right species resembles the query in the very "
                             "properties that vary within a genus.")
    parser.add_argument("--taxon-override", action="append", metavar="GENOME=TAXON",
                        help="Search TAXON instead of the MEGAN genus call for GENOME, e.g. "
                             "P2_G9=Enterobacteriaceae. Use when the genus call is a strain "
                             "identifier or an obsolete name with no reference genomes. "
                             "Accepts genus, family or order.")
    return parser.parse_args()


def main() -> None:
    options = parse_args()
    if options.genome:
        samples = list(options.genome)
    else:
        with options.manifest.open(encoding="utf-8") as handle:
            samples = [row["genome"] for row in csv.DictReader(handle, delimiter="\t")]

    overrides: dict[str, str] = {}
    for item in options.taxon_override or []:
        genome_name, _, taxon = item.partition("=")
        if not taxon:
            raise SystemExit(f"--taxon-override needs GENOME=TAXON, got {item!r}")
        overrides[genome_name] = taxon

    queries: list[tuple[str, Path, str, str, dict[str, str]]] = []
    for sample in samples:
        strain_dir = options.genome_root / sample
        if not strain_dir.is_dir() and (LOCAL_GENOME_ROOT / sample).is_dir():
            strain_dir = LOCAL_GENOME_ROOT / sample
            log(f"  {sample}: drive path unavailable, using the local panel copy")
        if not strain_dir.is_dir():
            log(f"WARNING: {sample}: folder not found under {options.genome_root} "
                f"or {LOCAL_GENOME_ROOT}")
            continue
        assembly = select_assembly(strain_dir)
        if assembly is None:
            log(f"WARNING: {sample}: no assembly FASTA found")
            continue
        genus, support = megan_genus(strain_dir)
        if sample in overrides:
            log(f"  {sample}: searching {overrides[sample]} instead of MEGAN call {genus!r}")
            genus, support = overrides[sample], support
        queries.append((sample, assembly, genus, support, quast_metrics(strain_dir),
                        megan_family(strain_dir)))

    genera = sorted({genus for _, _, genus, _, _, _ in queries if genus})
    log(f"{len(queries)} assemblies across genera: {', '.join(genera) or 'none'}")

    overridden = set(overrides.values())
    genus_tax_ids: dict[str, int] = {}
    for genus in genera:
        # An overridden taxon may be a family or order, so allow those ranks for it only.
        ranks = ("GENUS", "FAMILY", "ORDER") if genus in overridden else ("GENUS",)
        tax_id = resolve_taxon_tax_id(genus, ranks)
        if tax_id is None:
            log(f"  {genus}: could not resolve to a bacterial genus tax id; skipping")
            continue
        genus_tax_ids[genus] = tax_id
        log(f"  {genus}: tax id {tax_id}")

    # Size/GC of the queries in each genus, used to prefilter the reference set.
    genus_query_stats: dict[str, tuple[int, float]] = {}
    for sample, assembly, genus, support, quast, _ in queries:
        try:
            length = int(quast.get("assembly_bp") or 0)
            gc = float(quast.get("assembly_gc_percent") or 0)
        except ValueError:
            length, gc = 0, 0.0
        if genus and length and gc:
            genus_query_stats.setdefault(genus, (length, gc))

    # No reference sets are built up front any more. Searching per sample, in batches,
    # is what keeps peak memory to one batch instead of every genus at once -- and it is
    # why the cap and the size/GC window are no longer needed to make this run.
    families: dict[str, tuple[int, str] | None] = {}
    identity_rows: list[dict[str, object]] = []
    all_rows: list[dict[str, object]] = []

    def record_hits(sample: str, scope: str, hits: list[tuple[float, float, str, dict]]) -> None:
        for ani, fraction, accession, reference in hits:
            all_rows.append({
                "genome": sample,
                "search_scope": scope,
                "reference": accession,
                "species": reference.get("species", ""),
                "organism": reference.get("organism", ""),
                "ani_percent": round(ani, 4),
                "aligned_fraction": round(fraction, 4),
            })

    for sample, assembly, genus, support, quast, megan_fam in queries:
        row: dict[str, object] = {
            "genome": sample,
            "megan_genus": genus,
            "megan_genus_percent": support,
            "search_scope": "",
            "assembly_path": str(assembly),
            **quast,
            "best_species": "",
            "best_reference": "",
            "ani_percent": "",
            "aligned_fraction": "",
            "second_species": "",
            "second_ani_percent": "",
            "call": "no_reference_set",
        }
        hits: list[tuple[float, float, str, dict]] = []
        scope = ""

        # A genome with no MEGAN taxonomy at all never reaches the genus or family
        # branches below, because both are keyed on the genus call -- it simply falls
        # through as no_reference_set. That is what happened to the seven worm-collection
        # isolates, which have no diamond_megan directory: they occupied four of the
        # thirteen gene-content clusters and stayed nameless. --fallback-family gives them
        # a reference set to be searched against, chosen by the caller rather than guessed.
        # Fires for a missing genus *and* for a genus that will not resolve to a tax id.
        # The second case is not an edge case: `P1_A10` reads Raoultella, which was merged
        # into Klebsiella in 2021 and no longer exists at genus rank, so the lookup fails
        # and without this the strain is skipped entirely rather than searched.
        if (not genus or genus not in genus_tax_ids) and options.fallback_family:
            # ORDER is accepted as well as FAMILY, and it is usually the right scope. The
            # 2016 reclassification split the old broad Enterobacteriaceae into
            # Yersiniaceae, Morganellaceae, Hafniaceae, Pectobacteriaceae, Budviciaceae and
            # Erwiniaceae, so a family search misses fermentative enterobacteria that now
            # sit outside it. Five of the seven unassigned isolates best-hit
            # Enterobacteriaceae at 79-84% ANI with alignment as low as 0.24, which is not
            # a weak match but a statement that the reference set is the wrong one. This is
            # the fourth post-2016 reclassification to bite this project, after
            # Brucella/Ochrobactrum, Mammaliicoccus/Staphylococcus and Methylorubrum.
            fallback_tax_id = resolve_taxon_tax_id(options.fallback_family,
                                                   ("FAMILY", "ORDER"))
            if fallback_tax_id:
                log(f"  {sample}: no MEGAN taxonomy; searching fallback family "
                    f"{options.fallback_family} ({fallback_tax_id})")
                hits = search_references(options.fallback_family, fallback_tax_id,
                                         assembly, options)
                scope = f"fallback-family:{options.fallback_family}"
                record_hits(sample, scope, hits)

        if genus in genus_tax_ids:
            length, gc = genus_query_stats.get(genus, (0, 0.0))
            hits = search_references(genus, genus_tax_ids[genus], assembly, options, length, gc)
            scope = f"genus:{genus}"
            record_hits(sample, scope, hits)

        # Escalate to family when the genus set holds no species-level match. A best
        # hit at 80% ANI means the genus reference set does not contain this organism,
        # so the MEGAN genus call itself is the thing in doubt.
        if (not hits or hits[0][0] < SPECIES_ANI) and genus in genus_tax_ids:
            if genus not in families:
                families[genus] = resolve_family_tax_id(genus_tax_ids[genus])
            candidates: list[tuple[int, str]] = []
            family = families.get(genus)
            if family:
                candidates.append(family)
            # MEGAN's own family row, when it disagrees with the genus's family, is the
            # better bet: the genus row is a single-name row normalised to 100%, whereas
            # the family row is a real competition between names. P1_B3 reads
            # "Aneurinibacillus 100%" at genus but "Bacillaceae 99.98% vs
            # Aneurinibacillaceae 0.02%" at family, and ANI agrees with the family row.
            if megan_fam:
                megan_family_name, megan_family_percent = megan_fam
                if not family or megan_family_name.casefold() != family[1].casefold():
                    megan_family_tax_id = resolve_taxon_tax_id(megan_family_name, ("FAMILY",))
                    if megan_family_tax_id:
                        log(f"  {sample}: MEGAN's family row says {megan_family_name} "
                            f"({megan_family_percent:.2f}%), which disagrees with the family of "
                            f"genus {genus}; trying it first")
                        candidates.insert(0, (megan_family_tax_id, megan_family_name))
            for family_tax_id, family_name in candidates:
                log(f"  {sample}: no species match in {genus}; escalating to family "
                    f"{family_name} ({family_tax_id})")
                family_hits = search_references(
                    family_name or f"family_{family_tax_id}", family_tax_id, assembly, options
                )
                record_hits(sample, "family", family_hits)
                if family_hits and (not hits or family_hits[0][0] > hits[0][0]):
                    hits = family_hits
                    scope = f"family:{family_name}"
                # Stop as soon as one family yields a species-level match; otherwise try
                # the next candidate. Ordering matters, not just the break: adopting the
                # hits before breaking is what keeps a successful escalation from being
                # discarded on its way out of the loop.
                if family_hits and family_hits[0][0] >= SPECIES_ANI:
                    break

        if hits:
            ani, fraction, accession, reference = hits[0]
            row.update({
                "search_scope": scope,
                "best_species": reference.get("species", ""),
                "best_reference": f"{accession} ({reference.get('organism', '')})",
                "ani_percent": round(ani, 3),
                "aligned_fraction": round(fraction, 4),
                "call": classify(ani, fraction),
            })
            if len(hits) > 1:
                row["second_species"] = hits[1][3].get("species", "")
                row["second_ani_percent"] = round(hits[1][0], 3)
        elif genus in genus_tax_ids or scope.startswith("fallback-family:"):
            # A searched-but-unmatched genome must not report "no_reference_set", which
            # would read as "never looked" rather than "looked and found nothing".
            row["call"] = "no_ani_hit"
        identity_rows.append(row)

    options.output_dir.mkdir(parents=True, exist_ok=True)
    merge_rows(options.output_dir / "strain_identity.tsv", identity_rows,
               key=lambda r: str(r["genome"]))
    if all_rows:
        merge_rows(options.output_dir / "ani_all_hits.tsv", all_rows,
                   key=lambda r: f"{r['genome']}|{r['reference']}",
                   replace_scope=lambda r: str(r["genome"]),
                   sort_key=lambda r: (r["genome"], -float(r["ani_percent"])))

    print()
    for row in identity_rows:
        print(f"{row['genome']:8} {str(row['megan_genus']):16} -> {str(row['best_species']):34} "
              f"ANI {row['ani_percent'] or 'n/a':>8}  aln {row['aligned_fraction'] or 'n/a':>7}  "
              f"[{row['call']}] via {row['search_scope'] or 'n/a'}")
        if row["second_species"]:
            print(f"{'':8} next best: {row['second_species']} at {row['second_ani_percent']}%")
    print(f"\nwrote {options.output_dir / 'strain_identity.tsv'}")


if __name__ == "__main__":
    main()
