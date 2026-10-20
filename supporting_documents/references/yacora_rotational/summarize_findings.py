"""Freeze concise bibliographic identities and explicit access limitations."""
import hashlib
import html
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(name):
    return json.loads((ROOT/name).read_text())


def normal(text):
    return re.sub(r'\s+', '', html.unescape(re.sub('<[^>]*>', '', text))).replace('’', "'")


def main():
    findings = load('findings.json')
    selected = [row for row in findings['papers'] if row['year'] in (2023, 2025, 2026)
                and row['role'] != 'related_vibrational_CR_model_for_multiple_emission_bands']
    assert len(selected) == 3
    crossrefs = {'10.1016/j.jqsrt.2025.109414': 'crossref-fulcher-2025.json',
                '10.1088/1361-6595/ae6df8': 'crossref-guide-2026.json',
                '10.1140/epjd/s10053-023-00706-5': 'crossref-data-needs-2023.json'}
    page = load('official-ipp-publications.json')
    assert page['status'] == 200
    checks = []
    for row in selected:
        response = load(crossrefs[row['DOI']])
        metadata = response['JSON_metadata']
        assert response['status'] == 200 and metadata['DOI'] == row['DOI']
        assert normal(metadata['title'][0]) == normal(row['title'])
        assert metadata['published']['date-parts'][0][0] == row['year']
        assert [name['family'] for name in metadata['author']] == [name.split()[-1] for name in row['authors']]
        assert 'https://dx.doi.org/'+row['DOI'] in page['links']
        checks.append({'DOI': row['DOI'], 'title_authors_year_matched_registry': True,
                       'official_IPP_DOI_link_matched': True, 'registry_response_record': crossrefs[row['DOI']]})
    summary = {'schema': 'powerflow-yacora-rotational-literature-summary-v1',
        'status': 'verified_three_primary_citations', 'checked_utc_date': '2026-10-20',
        'primary_citations': selected, 'official_help': findings['official_help'],
        'official_author_publication_page': findings['official_author_publication_page'],
        'recommended_comparison': 'Name the exact Yacora model variant and data revision. Compare the2025H2Fulcher corona extension only within its band/process domain, separately from electronic-only Yacora on the Web and from the independently supplied plasma implementation.',
        'VUV_limit': 'Fulcher-α rotational resolution does not establish a rotationally resolved Lyman/Werner115–165nm observation model. The2023overview distinguishes selected UV corona models; use their actual state/data basis if required.',
        'coverage_limit': 'No H2O or complete molecular neutral/ion/anion network capability is established.',
        'access_limits': ['Current ordinary Yacora origin GETs returned503; official H2help was read through its public indexed content.',
            'Official IPP author publication page and all three exact DOI registry records returned200.',
            'Publisher-indexed abstract supports the2025Fulcher model scope; direct publisher page read returned403.',
            'No executable model, account calculation, data download, full-paper PDF, source code or runtime comparison was performed.'],
        'public_solver_data_or_document_edits': False, 'external_messages_sent': False}
    (ROOT/'summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False)+'\n')
    paths = ['summary.json','findings.json','requests.jsonl','indexed-query-records.json','official-ipp-publications.json',
             *crossrefs.values(), 'summarize_findings.py','read_official_pages.py',
             'official-home.json','official-home-apex.json','official-home-http.json']
    verification = {'schema': 'powerflow-yacora-literature-verification-v1', 'status': 'passed',
        'primary_citations_checked': 3, 'exact_registry_checks': checks,
        'artifact_sha256': {name:sha(ROOT/name) for name in paths},
        'verification_scope': 'Saved bibliographic metadata and source-association checks only; no scientific model validation.',
        'direct_official_Yacora_origin_currently_accessible': False,
        'full_copyrighted_PDFs_downloaded_or_committed': False,
        'CR_solver_data_or_public_document_edits': False, 'external_messages_sent': False}
    (ROOT/'verification.json').write_text(json.dumps(verification, indent=2, ensure_ascii=False)+'\n')
    print(json.dumps({'summary_sha256':sha(ROOT/'summary.json'),
                      'verification_sha256':sha(ROOT/'verification.json'), 'primary_citations':3},indent=2))


if __name__ == '__main__':
    main()
