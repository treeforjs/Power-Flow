"""Read two newly acquired catalog records without guessing document URLs."""
from html.parser import HTMLParser
import json
from pathlib import Path
from urllib.parse import urljoin

class Page(HTMLParser):
    def __init__(self):
        super().__init__()
        self.meta, self.links, self.text = [], [], []
    def handle_starttag(self, tag, attrs):
        values = dict(attrs)
        if tag == 'meta':
            self.meta.append(values)
        if tag == 'a' and values.get('href'):
            self.links.append(values['href'])
    def handle_data(self, data):
        self.text.append(data)

work = Path(__file__).resolve().parent
cache = Path('/workspace/research_sources/npl_processability_1996')
result = {}
for label, base in [('density_report_landing', 'https://eprintspublications.npl.co.uk/1034/'),
                    ('programme_paper_landing', 'https://eprintspublications.npl.co.uk/1324/')]:
    parser = Page()
    parser.feed((cache / (label + '.body')).read_text())
    result[label] = {
        'metadata': [item for item in parser.meta if item.get('name', '').lower().startswith(('citation', 'dc.'))],
        'document_links': [urljoin(base, link) for link in parser.links if '.pdf' in link.lower()],
        'visible_text': ' '.join(' '.join(parser.text).split()),
    }
out = work / 'landing-inspection.json'
assert not out.exists()
out.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
