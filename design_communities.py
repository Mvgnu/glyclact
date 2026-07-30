#!/usr/bin/env python3
"""Design the community x medium landscape for the glycerol -> lactate pool.

The feature vector is [p1..p3, s1..s8] with **at least one producer always present**,
crossed with medium conditions. The existing pipeline designer cannot express that: it
samples communities by size and either forces one focal species into all of them or
allows arbitrary subsets. Neither gives "any non-empty producer subset x any supporter
subset", which is the design here -- a community with no producer cannot make lactate, so
spending simulations on those measures nothing.

What the output is for: the **end-state composition**. Total biomass, the producers' share
of it and the supporters' share are what determine achievable lactate, because a producer
that gets outcompeted contributes nothing regardless of its monoculture capacity. That is
the quantity dFBA is determinate for -- unlike lactate flux, which is degenerate
(PRODUCT_PREDICTION_LIMITS.md).

Writes a manifest in the schema run_batch.py consumes.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
from pathlib import Path

import pandas as pd
import yaml

PRODUCERS = ["P1_B3", "P1_B11", "P3_F11"]
SUPPORTERS = ["P3_C11", "P3_B12", "P2_G9", "P3_F8", "P2_D5", "P2_D12", "P3_G8", "P2_E7"]


def stable_seed(*parts: object, modulo: int = 2_147_483_647) -> int:
    digest = hashlib.sha256("|".join(str(p) for p in parts).encode()).hexdigest()
    return int(digest, 16) % modulo


def flatten_medium(medium: dict[str, object]) -> dict[str, float]:
    flat: dict[str, float] = {}
    for group in medium.values():
        if isinstance(group, dict):
            for metabolite, concentration in group.items():
                flat[str(metabolite)] = float(concentration)
    return flat


def build_media(library_path: Path, wanted: list[str]) -> list[dict[str, object]]:
    library = yaml.safe_load(library_path.read_text(encoding="utf-8"))
    background = library.get("background_medium", {})
    definitions = library["media_definitions"]
    media: list[dict[str, object]] = []
    for name in wanted:
        if name not in definitions:
            raise SystemExit(f"medium {name!r} is not in {library_path}")
        combined = {**background, **definitions[name]}
        media.append({"medium": name, "definition": flatten_medium(combined)})
    return media


def communities(
    producers: list[str],
    supporters: list[str],
    max_supporters: int | None,
    min_producers: int,
    max_producers: int,
) -> list[tuple[str, ...]]:
    """Every non-empty producer subset crossed with every supporter subset.

    ``max_supporters`` truncates by community size rather than by sampling, so the design
    stays a complete factorial of whatever size range it covers -- which is what makes all
    the interaction terms estimable instead of aliased.
    """
    producer_sets = [
        combo
        for size in range(min_producers, max_producers + 1)
        for combo in itertools.combinations(producers, size)
    ]
    limit = len(supporters) if max_supporters is None else max_supporters
    supporter_sets = [
        combo
        for size in range(0, limit + 1)
        for combo in itertools.combinations(supporters, size)
    ]
    out = []
    for producer_set in producer_sets:
        for supporter_set in supporter_sets:
            out.append(tuple(sorted(producer_set + supporter_set)))
    return out


def main() -> None:
    root = Path(__file__).parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--media-library", type=Path, default=root / "comets" / "media_library.yaml")
    parser.add_argument("--output", type=Path, default=root / "comets" / "manifest.tsv")
    parser.add_argument("--media", nargs="+",
                        default=["glycerol_minimal", "glycerol_glucose",
                                 "glycerol_yeast_extract", "glycerol_glucose_yeast_extract"])
    parser.add_argument("--min-producers", type=int, default=1)
    parser.add_argument("--max-producers", type=int, default=3)
    parser.add_argument("--max-supporters", type=int, default=None,
                        help="Cap supporters per community. Omit for the full 2^8 landscape.")
    parser.add_argument("--initial-biomass", type=float, default=0.001)
    parser.add_argument("--replicates", type=int, default=1)
    options = parser.parse_args()

    media = build_media(options.media_library, options.media)
    groups = communities(PRODUCERS, SUPPORTERS, options.max_supporters,
                         options.min_producers, options.max_producers)

    rows = []
    index = 1
    for community in groups:
        present_producers = [s for s in community if s in PRODUCERS]
        for medium in media:
            for replicate in range(options.replicates):
                rows.append({
                    "sim_id": f"sim_{index:05d}",
                    "community": ",".join(community),
                    "community_size": len(community),
                    "n_producers": len(present_producers),
                    "n_supporters": len(community) - len(present_producers),
                    "medium": medium["medium"],
                    "environment_variable": "medium",
                    "environment_value": medium["medium"],
                    # run_batch keeps target_* diagnostics; point them at a present producer
                    # so they are well-defined. The labels used here are compositional.
                    "target_species": present_producers[0],
                    "initial_biomass": options.initial_biomass,
                    "replicate": replicate,
                    "random_seed": stable_seed("sim", index),
                    "medium_json": yaml.safe_dump(medium["definition"], sort_keys=True).strip(),
                })
                index += 1

    frame = pd.DataFrame(rows)
    options.output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(options.output, sep="\t", index=False)

    print(f"producers {PRODUCERS} (>=1 present), supporters {SUPPORTERS}")
    print(f"communities: {len(groups)}   media: {len(media)}   replicates: {options.replicates}")
    print(f"simulations: {len(rows)}")
    print(frame.groupby("community_size").size().to_string())
    print(f"\nwrote {options.output}")


if __name__ == "__main__":
    main()
