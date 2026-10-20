"""Resolve local publication links without opening network resources."""
from html import unescape
from html.parser import HTMLParser
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit


def prose(text, remove_inline_code=True):
    lines = []
    fence = None
    for line in text.splitlines(keepends=True):
        match = re.match(r'^\s{0,3}(`{3,}|~{3,})', line)
        if match:
            marker = match.group(1)
            if fence is None:
                fence = marker
            elif marker[0] == fence[0] and len(marker) >= len(fence):
                fence = None
            lines.append('\n' if line.endswith('\n') else '')
        elif fence is not None:
            lines.append('\n' if line.endswith('\n') else '')
        else:
            lines.append(line)
    cleaned = ''.join(lines)
    return re.sub(r'(`+)([^`]*?)\1', lambda m: ' '*len(m.group(0)), cleaned) if remove_inline_code else cleaned


def destination(raw):
    raw = raw.strip()
    if raw.startswith('<'):
        stop = raw.find('>')
        if stop < 0:
            raise ValueError('Unclosed angle-bracket link destination')
        raw = raw[1:stop]
    else:
        raw = re.split(r'\s', raw, maxsplit=1)[0]
    return unescape(re.sub(r'\\([\\`*{}\[\]()#+\-.!_> ])', r'\1', raw))


def _close(text, start, opening, closing):
    depth = 1
    i = start+1
    while i < len(text):
        if text[i] == '\\':
            i += 2
            continue
        if text[i] == opening:
            depth += 1
        elif text[i] == closing:
            depth -= 1
            if depth == 0:
                return i
        i += 1
    return None


def targets(text):
    cleaned = prose(text)
    references = {}
    scanned = []
    links = []
    lines = cleaned.splitlines(keepends=True)
    skipped = set()
    for line_number, line in enumerate(lines, 1):
        if line_number in skipped:
            scanned.append('\n' if line.endswith('\n') else '')
            continue
        match = re.match(r'^\s{0,3}\[([^\]]+)\]:\s*(<[^>]*>|\S+)', line)
        if match is None and re.match(r'^\s{0,3}\[[^\]]+\]:\s*$', line) and line_number < len(lines):
            match = re.match(r'^\s{0,3}\[([^\]]+)\]:\s*(<[^>]*>|\S+)', line.rstrip()+' '+lines[line_number].strip())
            if match:
                skipped.add(line_number+1)
        if match:
            label, target = match.groups()
            key = ' '.join(label.split()).casefold()
            references[key] = destination(target)
            links.append({'line': line_number, 'target': references[key], 'kind': 'reference_definition'})
            scanned.append('\n' if line.endswith('\n') else '')
        else:
            scanned.append(line)
    scanned = ''.join(scanned)
    i = 0
    while i < len(scanned):
        if scanned[i] != '[' or (i and scanned[i-1] == '\\'):
            i += 1
            continue
        end = _close(scanned, i, '[', ']')
        if end is None:
            i += 1
            continue
        label = scanned[i+1:end]
        after = end+1
        line_number = scanned.count('\n', 0, i)+1
        if '[' in label:
            links.extend({**nested, 'line': line_number+nested['line']-1} for nested in targets(label))
        if after < len(scanned) and scanned[after] == '(':
            stop = _close(scanned, after, '(', ')')
            if stop is None:
                raise ValueError(f'Unclosed inline link at line {line_number}')
            links.append({'line': line_number, 'target': destination(scanned[after+1:stop]), 'kind': 'inline'})
            i = stop+1
        elif after < len(scanned) and scanned[after] == '[':
            stop = _close(scanned, after, '[', ']')
            if stop is None:
                i = end+1
                continue
            key = ' '.join((scanned[after+1:stop] or label).split()).casefold()
            if key in references:
                links.append({'line': line_number, 'target': references[key], 'kind': 'reference'})
            i = stop+1
        else:
            key = ' '.join(label.split()).casefold()
            if key in references:
                links.append({'line': line_number, 'target': references[key], 'kind': 'shortcut_reference'})
            i = end+1
    return links


class _HTMLTargets(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []
        self.ids = set()

    def handle_starttag(self, tag, attributes):
        for key, value in attributes:
            if value is None:
                continue
            if key in {'href', 'src'}:
                self.links.append({'line': self.getpos()[0], 'target': value, 'kind': 'html_'+key})
            if key == 'id' or (tag == 'a' and key == 'name'):
                self.ids.add(value)


def anchors(text):
    cleaned = prose(text, remove_inline_code=False)
    parser = _HTMLTargets()
    parser.feed(prose(text))
    result, seen = set(parser.ids), {}
    lines = cleaned.splitlines()
    for i, line in enumerate(lines):
        match = re.match(r'^\s{0,3}#{1,6}\s+(.+?)(?:\s+#+\s*)?$', line)
        heading = match.group(1) if match else None
        if heading is None and i+1 < len(lines) and line.strip() and re.fullmatch(r'\s{0,3}(?:=+|-+)\s*', lines[i+1]):
            heading = line.strip()
        if heading is None:
            continue
        heading = re.sub(r'!?\[([^\]]+)\]\([^)]*\)', r'\1', heading)
        heading = re.sub(r'(`+)(.*?)\1', r'\2', heading)
        heading = re.sub(r'(?<!\w)(\*{1,3}|_{1,3})(.+?)\1(?!\w)', r'\2', heading)
        heading = re.sub(r'<[^>]*>', '', unescape(heading)).lower()
        slug = re.sub(r'[^\w\- ]', '', heading).replace(' ', '-')
        counter = seen.get(slug, 0)
        seen[slug] = counter+1
        result.add(slug if counter == 0 else slug+'-'+str(counter))
    return result


def verify_documents(repo, documents, *, allowed_pending_outputs=()):
    repo = Path(repo).resolve()
    allowed_pending_outputs = {Path(path).resolve() for path in allowed_pending_outputs}
    result, anchor_cache = [], {}
    for document in documents:
        text = document.read_text()
        html = _HTMLTargets()
        html.feed(prose(text))
        for item in targets(text)+html.links:
            uri = urlsplit(item['target'])
            if uri.scheme or uri.netloc:
                if uri.scheme.lower() == 'file':
                    raise ValueError(f'Local file URI in {document}: {item}')
                continue
            path = unquote(uri.path)
            if path.startswith('/') or '\x00' in path:
                raise ValueError(f'Nonportable local link in {document}: {item}')
            target = (document.parent/path).resolve() if path else document.resolve()
            pending = target in allowed_pending_outputs and not target.exists()
            if not target.is_relative_to(repo) or (not target.exists() and not pending):
                raise ValueError(f'Broken or escaping local link in {document}: {item}')
            fragment = unquote(uri.fragment)
            anchor_target = target
            if pending and fragment:
                raise ValueError(f'Unverifiable fragment on planned output in {document}: {item}')
            if fragment and target.is_dir():
                readmes = [p for p in target.iterdir() if p.is_file() and p.name.lower() == 'readme.md']
                if len(readmes) != 1:
                    raise ValueError(f'Directory fragment has no unique README in {document}: {item}')
                anchor_target = readmes[0]
            if fragment and anchor_target.suffix.lower() in {'.md', '.markdown'}:
                if anchor_target not in anchor_cache:
                    anchor_cache[anchor_target] = anchors(anchor_target.read_text())
                if fragment not in anchor_cache[anchor_target] and fragment.removeprefix('user-content-') not in anchor_cache[anchor_target]:
                    raise ValueError(f'Unknown Markdown anchor in {document}: {item}')
            elif fragment and re.fullmatch(r'L\d+(?:-L\d+)?', fragment):
                indices = [int(part.removeprefix('L')) for part in fragment.split('-')]
                if not 1 <= min(indices) <= max(indices) <= len(target.read_text().splitlines()):
                    raise ValueError(f'Unknown line anchor in {document}: {item}')
            result.append({'document': str(document.relative_to(repo)), **item,
                           'resolved_target': str(target.relative_to(repo)), 'fragment': fragment,
                           'planned_output': pending})
    return result
