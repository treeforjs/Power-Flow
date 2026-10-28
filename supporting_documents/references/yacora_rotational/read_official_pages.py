"""Record ordinary public-page requests and concise bibliographic evidence."""
import argparse
import datetime
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import urllib.error
import urllib.parse
import urllib.request

ROOT = Path(__file__).resolve().parent


class Page(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []
        self.text = []
        self.hidden = 0

    def handle_starttag(self, tag, attrs):
        if tag in ('script', 'style'):
            self.hidden += 1
        if tag == 'a':
            href = dict(attrs).get('href')
            if href:
                self.links.append(href)
        if tag in ('p', 'div', 'br', 'li', 'h1', 'h2', 'h3', 'tr'):
            self.text.append('\n')

    def handle_endtag(self, tag):
        if tag in ('script', 'style'):
            self.hidden -= 1

    def handle_data(self, data):
        if not self.hidden:
            self.text.append(data)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('url')
    parser.add_argument('--label', required=True)
    args = parser.parse_args()
    assert urllib.parse.urlsplit(args.url).scheme in ('http', 'https')
    stamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
    request = urllib.request.Request(args.url, headers={'User-Agent': 'Research bibliographic verification'})
    try:
        with urllib.request.urlopen(request, timeout=35) as response:
            body = response.read(4*1024*1024+1)
            status, headers, final_url = response.status, dict(response.headers), response.geturl()
        assert len(body) <= 4*1024*1024
        error = None
    except urllib.error.HTTPError as exc:
        body = exc.read(1024*1024)
        status, headers, final_url, error = exc.code, dict(exc.headers), exc.geturl(), str(exc)
    except Exception as exc:
        body, status, headers, final_url, error = b'', None, {}, args.url, str(exc)
    record = {'requested_utc': stamp, 'request_url': args.url, 'method': 'GET',
        'final_url': final_url, 'status': status, 'bytes': len(body),
        'sha256': hashlib.sha256(body).hexdigest(),
        'content_type': next((v for k,v in headers.items() if k.lower() == 'content-type'), ''), 'error': error,
        'proxy_and_CA_settings': 'Inherited unchanged; no bypass or TLS override.',
        'response_headers': {k:v for k,v in headers.items() if k.lower() in
                             ('content-type','content-length','last-modified','etag','location')}}
    with (ROOT/'requests.jsonl').open('a') as stream:
        stream.write(json.dumps(record)+'\n')
    page = Page()
    if 'html' in record['content_type']:
        page.feed(body.decode('utf-8', errors='replace'))
        text = '\n'.join(' '.join(line.split()) for line in ''.join(page.text).splitlines() if line.strip())
        evidence = {**record, 'visible_text': text[:22000],
            'visible_text_truncated': len(text)>22000,
            'links': sorted(set(urllib.parse.urljoin(final_url, x) for x in page.links))}
    elif 'json' in record['content_type']:
        document = json.loads(body)
        message = document.get('message', document)
        fields = ('DOI','title','author','published','published-print','published-online',
                  'container-title','volume','issue','page','article-number','URL','link')
        evidence = {**record, 'JSON_metadata': {key:message[key] for key in fields if key in message}}
    else:
        evidence = record
    path = ROOT/(args.label+'.json')
    assert not path.exists()
    path.write_text(json.dumps(evidence, indent=2, ensure_ascii=False)+'\n')
    print(json.dumps(evidence, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
