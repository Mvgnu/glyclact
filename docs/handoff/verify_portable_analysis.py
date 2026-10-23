"""Verify October analysis in a disposable copy; never rewrite frozen results."""
import ast
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
ECHO = ROOT/'protocols/echo'


def main():
    original = json.loads((ECHO/'readout_analysis/20261002_diagnostic/analysis.json').read_text())
    with tempfile.TemporaryDirectory(prefix='lactate-pro-portability-') as tmp:
        echo = (Path(tmp)/'protocols/echo').resolve()
        shutil.copytree(ECHO, echo, ignore=shutil.ignore_patterns('__pycache__', '.DS_Store', '._*'))
        subprocess.run([sys.executable, str(echo/'analyze_readout_20261002.py')], check=True, capture_output=True, text=True)
        generated_path = echo/'readout_analysis/20261002_diagnostic/analysis.json'
        generated = json.loads(generated_path.read_text())
        for key in ('condition_summary', 'culture_pairs', 'duplicate_agreement',
                    'read_comparison', 'standard_diagnostics', 'postshake_standard_values',
                    'spike_pairs', 'background_pairs'):
            assert original[key] == generated[key], key
        # Emulate a different computer with stale absolute provenance paths.
        generated['design_path'] = '/old-host/protocols/echo/' + str(
            Path(generated['design_path']).relative_to(echo))
        for source in generated['sources']:
            source['file'] = '/old-host/Downloads/read-02-10/' + Path(source['file']).name
        generated_path.write_text(json.dumps(generated))
        subprocess.run([sys.executable, str(echo/'visualize_readout_20261002.py')], check=True)
    # Frozen timestamped files must still match their existing manifest entries.
    count = 0
    for manifest_path in (ECHO/'experiment_setup_output').rglob('manifest.json'):
        manifest = json.loads(manifest_path.read_text())
        files = manifest.get('files', {})
        if not isinstance(files, dict):
            continue
        for name, expected in files.items():
            if not isinstance(expected, str) or len(expected) != 64:
                continue
            path = manifest_path.parent/name
            if not path.is_file():
                continue  # Other historic manifests can reference external inputs.
            assert hashlib.sha256(path.read_bytes()).hexdigest() == expected, str(path)
            count += 1
    nb = json.loads((ECHO/'experiment_setup.ipynb').read_text())
    for i, cell in enumerate(nb['cells']):
        if cell['cell_type'] == 'code':
            ast.parse(''.join(cell['source']), filename=f'notebook-cell-{i}')
    print(f'Portable numerical results unchanged; figure paths portable; {count} frozen hashes verified; notebook parses.')


if __name__ == '__main__':
    main()
