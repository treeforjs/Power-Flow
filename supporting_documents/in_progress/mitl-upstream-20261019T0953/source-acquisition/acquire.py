"""Acquire one bounded public repository response without executing its content."""
import datetime
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
from urllib.parse import urlsplit

label, url = sys.argv[1:]
assert re.fullmatch('[a-z0-9_-]{1,60}', label)
assert urlsplit(url).scheme == 'https'
assert urlsplit(url).hostname in {'eprints.npl.co.uk', 'eprintspublications.npl.co.uk', 'www.npl.co.uk'}
cache = Path('/workspace/research_sources/npl_processability_1996')
work = Path(__file__).resolve().parent
body, headers = cache / (label + '.body'), cache / (label + '.headers')
record = work / (label + '.json')
assert not any(path.exists() for path in (body, headers, record))
command = ['curl', '--location', '--max-redirs', '5', '--connect-timeout', '12',
           '--max-time', '45', '--max-filesize', '12582912', '--proto', '=https',
           '--proto-redir', '=https', '--silent', '--show-error', '--dump-header', str(headers),
           '--output', str(body), '--write-out', '%{json}', url]
started = datetime.datetime.now(datetime.timezone.utc).isoformat()
result = subprocess.run(command, capture_output=True, text=True)
try:
    response = json.loads(result.stdout)
except ValueError:
    response = {'unparsed_stdout': result.stdout}
files = {str(path): {'bytes': path.stat().st_size, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
         for path in (body, headers) if path.exists()}
entry = {'schema': 'mitl-bounded-source-response-v1', 'started_utc': started,
         'finished_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
         'command': command, 'exit_code': result.returncode, 'stderr': result.stderr,
         'response': response, 'files': files, 'body_is_PDF': body.exists() and body.read_bytes().startswith(b'%PDF-'),
         'parameter_adoption': False}
record.write_text(json.dumps(entry, indent=2) + '\n')
print(json.dumps({'record': str(record), 'exit_code': result.returncode,
                  'http_code': response.get('http_code'), 'effective_url': response.get('url_effective'),
                  'content_type': response.get('content_type'), 'files': files}, indent=2))
raise SystemExit(result.returncode)
