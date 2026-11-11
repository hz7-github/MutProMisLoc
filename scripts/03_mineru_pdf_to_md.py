"""Install MinerU 4.x, convert PDFs and remove embedded image base64 from Markdown.

Install only: python scripts/03_mineru_pdf_to_md.py --install
Preview: python scripts/03_mineru_pdf_to_md.py --plan-only
Convert: python scripts/03_mineru_pdf_to_md.py --limit 10
Re-clean existing raw Markdown: add --clean-only.
"""
from __future__ import annotations

import argparse
import base64
import contextlib
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import signal
import sqlite3
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
FIELDS = ['source', 'raw_md', 'clean_md', 'status', 'removed_images', 'clean_chars', 'error', 'updated_at']
CLEANER_VERSION = '2'


def stop_process_tree(process):
    """Stop the dedicated child tree before temporary outputs are removed."""
    previous = signal.signal(signal.SIGINT, signal.SIG_IGN)
    try:
        if os.name == 'nt':
            if process.poll() is None:
                result = subprocess.run(
                    ['taskkill', '/PID', str(process.pid), '/T', '/F'],
                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                    timeout=15, check=False)
                if result.returncode and process.poll() is None:
                    raise RuntimeError('Could not stop the MinerU process tree.')
        else:
            # Children share this private process group, even if its leader exits.
            with contextlib.suppress(ProcessLookupError):
                os.killpg(process.pid, signal.SIGKILL)
        process.wait(timeout=15)
    finally:
        signal.signal(signal.SIGINT, previous)


def run_managed_command(command, *, timeout=None, check=True, **kwargs):
    """Own the external process tree on interrupt and timeout."""
    isolation = ({'creationflags': subprocess.CREATE_NEW_PROCESS_GROUP}
                 if os.name == 'nt' else {'start_new_session': True})
    process = subprocess.Popen(command, **isolation, **kwargs)
    try:
        # An infinite Windows process wait can defer Python's SIGINT handler.
        # Return to the interpreter frequently, even with no overall timeout.
        deadline = time.monotonic() + timeout if timeout is not None else None
        while (code := process.poll()) is None:
            remaining = deadline - time.monotonic() if deadline is not None else None
            if remaining is not None and remaining <= 0:
                raise subprocess.TimeoutExpired(command, timeout)
            time.sleep(min(0.1, remaining) if remaining is not None else 0.1)
    except (KeyboardInterrupt, subprocess.TimeoutExpired):
        print('Stopping external process tree; please wait.', flush=True)
        stop_process_tree(process)
        raise
    if check and code:
        raise subprocess.CalledProcessError(code, command)
    return subprocess.CompletedProcess(command, code)


IMAGE_URI = re.compile(
    r'data:image/[\w.+-]+;base64,[A-Za-z0-9+/=]+'
    r'(?:\r?\n(?:[A-Za-z0-9+/=]{64,}|[A-Za-z0-9+/=]+(?=["\x27)>])))*'
    r'(?=$|[\s"\x27)>])', re.I)


def clean_with_audit(text):
    """Remove image nodes with stable IDs and retain source/clean character spans."""
    text = text.lstrip('\ufeff').replace('\r\n', '\n').replace('\r', '\n')
    nodes = []
    patterns = [
        ('html_image', r'<img\b[^>]*\bsrc=["\x27]data:image/[^"\x27]+["\x27][^>]*>'),
        ('inline_markdown', r'!\[(?:\\.|[^\]])*\]\(\s*data:image/[^)]*\)'),
    ]
    for kind, pattern in patterns:
        for match in re.finditer(pattern, text, re.I):
            if not IMAGE_URI.search(match.group()):
                continue
            if not any(match.start() < end and match.end() > start for start, end, _ in nodes):
                nodes.append((match.start(), match.end(), kind))
    for match in IMAGE_URI.finditer(text):
        if not any(match.start() < end and match.end() > start for start, end, _ in nodes):
            line_start = text.rfind('\n', 0, match.start()) + 1
            kind = 'reference_definition' if re.match(r'\s*\[[^\]]+\]:\s*$', text[line_start:match.start()]) else 'bare_data_uri'
            nodes.append((match.start(), match.end(), kind))
    parts, images, cursor, offset = [], [], 0, 0
    for start, end, kind in sorted(nodes):
        prefix = text[cursor:start]
        parts.append(prefix)
        offset += len(prefix)
        node = text[start:end]
        uri = IMAGE_URI.search(node).group()
        payload = re.sub(r'\s+', '', uri.split(',', 1)[1])
        alt = ''
        if kind == 'inline_markdown':
            alt = re.match(r'!\[((?:\\.|[^\]])*)\]', node).group(1)
        elif kind == 'html_image':
            match = re.search(r'\balt=(["\x27])(.*?)\1', node, re.I | re.S)
            alt = match.group(2) if match else ''
        image_id = f'img{len(images)+1:04d}'
        replacement = f'[IMAGE_REMOVED id={image_id}]' + (' ' + alt if alt else '') + '\n' * node.count('\n')
        decoded_hash, decode_error = '', ''
        try:
            decoded_hash = hashlib.sha256(base64.b64decode(payload, validate=True)).hexdigest()
        except ValueError:
            decode_error = 'Invalid base64 encoding'
        images.append(dict(image_id=image_id, kind=kind, alt=alt, original_offset=[start, end],
                           clean_offset=[offset, offset+len(replacement)], replacement=replacement,
                           mime=uri[5:].split(';', 1)[0], encoded_bytes=len(payload), sha256=decoded_hash,
                           decode_error=decode_error, review_status='UNREVIEWED'))
        parts.append(replacement)
        offset += len(replacement)
        cursor = end
    parts.append(text[cursor:])
    cleaned = ''.join(parts)
    if re.search(r'data:image/[^,\s]*;base64,', cleaned, re.I):
        raise ValueError('Unparsed image base64 remains; inspect the raw Markdown')
    if not cleaned.strip():
        raise ValueError('Markdown is empty after cleaning')
    units, section = [], 'Document'
    for match in re.finditer(r'\S(?:.*?\S)?(?=\n[ \t]*\n|\s*\Z)', cleaned, re.S):
        heading = re.match(r'#{1,6}\s+([^\n]+)', match.group())
        if heading:
            section = heading.group(1)
        units.append(dict(id=f'p{len(units)+1:04d}', section=section, text=match.group(),
                          offset=[match.start(), match.end()]))
    return cleaned, images, units


def clean_markdown(text):
    cleaned, images, _ = clean_with_audit(text)
    return cleaned, len(images)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def atomic_text(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', dir=path.parent,
                                     suffix='.tmp', delete=False) as stream:
        temp = Path(stream.name)
        stream.write(text)
    try:
        os.replace(temp, path)
    finally:
        temp.unlink(missing_ok=True)


def environment_python(directory):
    return directory / ('Scripts/python.exe' if os.name == 'nt' else 'bin/python')


def install_mineru(args):
    if not (3, 11) <= sys.version_info[:2] < (3, 15):
        raise ValueError('MinerU installation requires Python 3.11 through 3.14')
    python = environment_python(args.mineru_env)
    commands = []
    if not python.is_file():
        commands.append([sys.executable, '-m', 'venv', str(args.mineru_env)])
    commands.append([str(python), '-m', 'pip', 'install', 'mineru>=4.0,<5'])
    for command in commands:
        print(subprocess.list2cmdline(command), flush=True)
        if not args.plan_only:
            run_managed_command(command, check=True)
    if not args.plan_only:
        run_managed_command([str(python), '-c',
                        "from importlib.metadata import version; print('MinerU', version('mineru'))"], check=True)
        print('Installation complete. Run this script again without --install to convert PDFs.')


def mineru_command(args):
    if args.mineru_command:
        return args.mineru_command
    local = args.mineru_env / ('Scripts/mineru-kit.exe' if os.name == 'nt' else 'bin/mineru-kit')
    return str(local) if local.is_file() else shutil.which('mineru-kit') or 'mineru-kit'


def build_command(args, source, destination):
    executable = mineru_command(args)
    launcher = [executable]
    if os.name == 'nt':
        resolved = Path(shutil.which(executable) or executable)
        python = resolved.parent / 'python.exe'
        if python.is_file():
            # Conda-based venvs may load an older MSVC DLL from the base Python
            # directory. Preload the system runtime only inside this child process.
            bootstrap = (
                "import ctypes,os,sys; from pathlib import Path; "
                "dll=Path(os.environ.get('SystemRoot',r'C:\\Windows'))/'System32'/'msvcp140.dll'; "
                "runtime=ctypes.WinDLL(str(dll)) if dll.is_file() else None; "
                "from mineru.kit.main import main; main()"
            )
            launcher = [str(python), '-c', bootstrap]
    command = launcher + ['parse', str(source), '-o', str(destination), '--tier', args.tier]
    if args.ocr_mode:
        command.extend(['--ocr-mode', args.ocr_mode])
    return command


def build_environment(args):
    env = os.environ.copy()
    for option, variable in [('small_backend', 'MINERU_MODEL_SMALL_BACKEND'),
                             ('vlm_engine', 'MINERU_MODEL_VLM_ENGINE'),
                             ('model_source', 'MINERU_MODEL_SOURCE')]:
        if getattr(args, option):
            env[variable] = getattr(args, option)
    return env


@contextlib.contextmanager
def conversion_lock(out):
    out.mkdir(parents=True, exist_ok=True)
    with (out / '.conversion.lock').open('a+b') as handle:
        handle.seek(0)
        handle.write(b'0')
        handle.flush()
        handle.seek(0)
        if os.name == 'nt':
            import msvcrt
            acquire = lambda: msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            release = lambda: msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
        else:
            import fcntl
            acquire = lambda: fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            release = lambda: fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
        try:
            acquire()
        except OSError as exc:
            raise ValueError('Another conversion process is using this output directory') from exc
        try:
            yield
        finally:
            handle.seek(0)
            release()


def convert_one(args, source, raw, clean, log, source_hash):
    if args.clean_only:
        content = source.read_text(encoding='utf-8-sig')
    else:
        raw.parent.mkdir(parents=True, exist_ok=True)
        log.parent.mkdir(parents=True, exist_ok=True)
        error = ''
        with tempfile.TemporaryDirectory(prefix='.mineru_', dir=raw.parent) as temporary:
            staged = Path(temporary) / raw.name
            for attempt in range(args.retries + 1):
                staged.unlink(missing_ok=True)
                try:
                    with log.open('a', encoding='utf-8') as stream:
                        stream.write(f'\nAttempt {attempt + 1}\n')
                        stream.flush()
                        run_managed_command(build_command(args, source, staged), stdout=stream,
                                       stderr=subprocess.STDOUT, env=build_environment(args),
                                       timeout=args.timeout or None, check=True)
                    content = staged.read_text(encoding='utf-8-sig')
                    if not content.strip():
                        raise ValueError('MinerU produced empty Markdown')
                    break
                except (subprocess.SubprocessError, OSError, ValueError) as exc:
                    error = f'{type(exc).__name__}: {exc}'
            else:
                raise ValueError(error)
            if digest(source) != source_hash:
                raise ValueError('Source changed during conversion; result not committed')
            # Keep the original text even if subsequent payload cleaning needs review.
            atomic_text(raw, content)
    cleaned, images, units = clean_with_audit(content)
    if digest(source) != source_hash:
        raise ValueError('Source changed during cleaning; result not committed')
    atomic_text(clean, cleaned)
    clean_hash = digest(clean)
    atomic_text(clean.with_suffix('.md.image_manifest.json'), json.dumps({
        'version': CLEANER_VERSION, 'source_file': str(raw.resolve()),
        'source_sha256': digest(raw), 'clean_sha256': clean_hash,
        'offset_convention': 'Half-open Unicode offsets after BOM removal and LF normalization',
        'images': images, 'statistics': {'images_found': len(images), 'raw_chars': len(content),
                                      'clean_chars': len(cleaned)}}, ensure_ascii=False, indent=2))
    atomic_text(clean.with_suffix('.md.source_units.json'), json.dumps({
        'clean_sha256': clean_hash, 'units': units}, ensure_ascii=False, indent=2))
    return len(images), len(cleaned)


def export_status(out, rows):
    import io
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=FIELDS, extrasaction='ignore')
    writer.writeheader()
    writer.writerows(rows)
    try:
        atomic_text(out / 'conversion_status.csv', '\ufeff' + stream.getvalue())
    except PermissionError:
        print('Close conversion_status.csv; progress remains saved in SQLite.', flush=True)


def run(args):
    if args.install:
        install_mineru(args)
        return 0
    source_root = args.raw_dir if args.clean_only else args.input_dir
    if not source_root.is_dir():
        raise ValueError('Input directory does not exist: ' + str(source_root))
    suffix = '.md' if args.clean_only else '.pdf'
    files = sorted((p for p in source_root.rglob('*') if p.is_file() and p.suffix.lower() == suffix),
                   key=lambda p: str(p).casefold())
    keys = [p.relative_to(source_root).with_suffix('.md').as_posix().casefold() for p in files]
    if len(keys) != len(set(keys)):
        raise ValueError('Input names would overwrite the same Markdown output')
    if args.plan_only:
        for source in files[:args.limit] if args.limit else files:
            raw = args.raw_dir / source.relative_to(source_root).with_suffix('.md')
            print('CLEAN ' + str(source) if args.clean_only else subprocess.list2cmdline(build_command(args, source, raw)))
        print('Clean Markdown:', args.output_dir / 'clean')
        return 0
    if not files:
        print('No matching input files; no conversion performed.')
        return 0
    if not args.clean_only and not shutil.which(mineru_command(args)):
        raise ValueError('MinerU CLI not found. Run this script with --install, or supply --mineru-command.')
    config = {key: getattr(args, key) for key in ('tier', 'ocr_mode', 'small_backend', 'vlm_engine', 'model_source')}
    config.update(cleaner=CLEANER_VERSION, mode='clean' if args.clean_only else 'convert', command=mineru_command(args))
    config_hash = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()
    rows = []
    attempts = 0
    failures = 0
    with conversion_lock(args.output_dir), contextlib.closing(sqlite3.connect(args.output_dir / '.conversion_state.sqlite3')) as db:
        db.execute('CREATE TABLE IF NOT EXISTS records (source TEXT PRIMARY KEY, result TEXT)')
        db.commit()
        for source in files:
            relative = source.relative_to(source_root).with_suffix('.md')
            raw = args.raw_dir / relative
            clean = args.output_dir / 'clean' / relative
            if source.resolve() == clean.resolve() or raw.resolve() == clean.resolve():
                raise ValueError('Raw/source and cleaned Markdown must use different paths')
            source_hash = digest(source)
            saved = db.execute('SELECT result FROM records WHERE source=?', (str(source.resolve()),)).fetchone()
            old = json.loads(saved[0]) if saved else {}
            row = dict(source=str(source.resolve()), raw_md=str(raw.resolve()), clean_md=str(clean.resolve()),
                       status='PENDING', error='', source_hash=source_hash, config_hash=config_hash)
            valid = (old.get('source_hash') == source_hash and old.get('config_hash') == config_hash
                     and old.get('status') in ('SUCCESS', 'SKIPPED') and raw.is_file() and clean.is_file()
                     and digest(raw) == old.get('raw_hash') and digest(clean) == old.get('clean_hash'))
            for suffix, field in (('.md.image_manifest.json', 'image_manifest_hash'),
                                  ('.md.source_units.json', 'source_units_hash')):
                sidecar = clean.with_suffix(suffix)
                valid = valid and sidecar.is_file() and digest(sidecar) == old.get(field)
            if valid and not args.force:
                row = dict(old, status='SKIPPED')
            elif args.limit is not None and attempts >= args.limit:
                rows.append(row)
                continue
            else:
                attempts += 1
                # A failed rerun must not leave an obsolete clean file for Step 04.
                clean.unlink(missing_ok=True)
                for suffix in ('.md.image_manifest.json', '.md.source_units.json'):
                    clean.with_suffix(suffix).unlink(missing_ok=True)
                try:
                    count, chars = convert_one(args, source, raw, clean,
                                                args.output_dir / 'logs' / relative.with_suffix('.log'), source_hash)
                    row.update(status='SUCCESS', removed_images=count, clean_chars=chars,
                               raw_hash=digest(raw), clean_hash=digest(clean),
                               image_manifest_hash=digest(clean.with_suffix('.md.image_manifest.json')),
                               source_units_hash=digest(clean.with_suffix('.md.source_units.json')))
                except (OSError, ValueError) as exc:
                    row.update(status='FAILED', error=str(exc))
                    failures += 1
                print(f"[{attempts}] {relative}: {row['status']}", flush=True)
            row['updated_at'] = datetime.now(timezone.utc).isoformat()
            db.execute('INSERT OR REPLACE INTO records VALUES (?, ?)', (str(source.resolve()), json.dumps(row)))
            db.commit()
            rows.append(row)
            export_status(args.output_dir, rows)
        export_status(args.output_dir, rows)
    print(f'Finished: {attempts} attempted; {failures} failed. Clean Markdown: {args.output_dir / "clean"}')
    return 1 if failures else 0


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--preset', choices=('original', 'expanded'), default='expanded')
    p.add_argument('--input-dir', type=Path)
    p.add_argument('--output-dir', type=Path)
    p.add_argument('--raw-dir', type=Path, help='Existing raw Markdown for --clean-only; default: <output-dir>/raw')
    p.add_argument('--install', action='store_true', help='Install MinerU in an isolated environment, then exit')
    p.add_argument('--mineru-env', type=Path, default=ROOT / '.mineru')
    p.add_argument('--mineru-command', help='Explicit path to an existing mineru-kit CLI')
    p.add_argument('--tier', choices=('flash', 'basic', 'standard', 'advanced'), default='standard')
    p.add_argument('--ocr-mode', choices=('txt', 'ocr'))
    p.add_argument('--small-backend', choices=('auto', 'onnx', 'torch'))
    p.add_argument('--vlm-engine', choices=('auto', 'llama-cpp', 'vllm', 'lmdeploy'))
    p.add_argument('--model-source', choices=('auto', 'huggingface', 'modelscope', 'local'))
    p.add_argument('--timeout', type=int, default=0, help='Per-paper timeout in seconds; 0 disables it')
    p.add_argument('--retries', type=int, default=1)
    p.add_argument('--limit', type=int, help='Maximum pending attempts; existing verified results are skipped')
    p.add_argument('--force', action='store_true', help='Regenerate existing conversion/cleaning results')
    p.add_argument('--clean-only', action='store_true', help='Clean raw Markdown without invoking MinerU')
    p.add_argument('--plan-only', '--dry-run', dest='plan_only', action='store_true')
    args = p.parse_args(argv)
    if args.timeout < 0 or args.retries < 0 or (args.limit is not None and args.limit < 1):
        p.error('Timeout/retries must be nonnegative; limit must be positive')
    args.input_dir = (args.input_dir or ROOT / 'data/02_download_papers' / args.preset).resolve()
    args.output_dir = (args.output_dir or ROOT / 'data/03_mineru_pdf_to_md' / args.preset).resolve()
    args.raw_dir = (args.raw_dir or args.output_dir / 'raw').resolve()
    return args


if __name__ == '__main__':
    try:
        raise SystemExit(run(parse_args()))
    except KeyboardInterrupt:
        print('Interrupted; completed conversions retained. Repeat the command to resume.', file=sys.stderr)
        raise SystemExit(130)
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        print(f'Error: {exc}', file=sys.stderr)
        raise SystemExit(2)
