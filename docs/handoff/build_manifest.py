"""Hash Git-eligible project files without following ignored dependency trees.

Run from any working directory. The manifest excludes itself, Git internals and
ignored runtime/cache/credential files. It records bytes, not inferred data quality.
"""
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT/'docs/handoff/PROJECT_FILE_MANIFEST.json'


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    listing = subprocess.check_output(
        ['git', 'ls-files', '--cached', '--others', '--exclude-standard', '-z'], cwd=ROOT)
    records = []
    for relative in sorted(set(x.decode() for x in listing.split(b'\0') if x)):
        path = ROOT/relative
        if path == OUTPUT or not path.exists():
            continue
        if path.is_symlink():
            raise ValueError(f'Review symlink before packaging: {relative}')
        records.append(dict(path=relative, bytes=path.stat().st_size, sha256=digest(path)))
    output = dict(snapshot_date='2026-10-23',
        scope='Available Git-eligible files in the glycerol_lactate_community repository',
        excluded='This manifest itself; Git internals; ignored dependencies, caches, credentials and OS metadata',
        file_count=len(records), total_bytes=sum(r['bytes'] for r in records), files=records)
    OUTPUT.write_text(json.dumps(output, indent=2)+'\n')
    print(f'{len(records)} files, {output["total_bytes"]:,} bytes: {OUTPUT}')


if __name__ == '__main__':
    main()
