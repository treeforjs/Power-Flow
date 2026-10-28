"""Extract candidate repository links from the saved official catalog response."""
from html.parser import HTMLParser
import json
from pathlib import Path
from urllib.parse import urljoin

class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.current = None
        self.rows = []
    def handle_starttag(self, tag, attrs):
        if tag == 'a':
            self.current = [dict(attrs).get('href', ''), []]
    def handle_data(self, data):
        if self.current is not None:
            self.current[1].append(data)
    def handle_endtag(self, tag):
        if tag == 'a' and self.current is not None:
            href, chunks = self.current
            self.rows.append({'url': urljoin('https://eprints.npl.co.uk/', href),
                              'text': ' '.join(''.join(chunks).split())})
            self.current = None

source = Path('/workspace/research_sources/npl_processability_1996/materials_catalog.body')
parser = Links()
parser.feed(source.read_text())
terms = ('levitated drop', 'processability', 'properties of alloys relevant',
         'solidification of commercial', 'modelling of heat and fluid')
rows = [row for row in parser.rows if any(term in row['text'].lower() for term in terms)]
out = Path(__file__).with_name('catalog-links.json')
assert not out.exists()
out.write_text(json.dumps({'source': str(source), 'selection_terms': terms, 'matches': rows}, indent=2) + '\n')
print(json.dumps(rows, indent=2))
