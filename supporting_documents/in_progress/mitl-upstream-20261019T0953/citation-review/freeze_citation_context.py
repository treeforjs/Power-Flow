"""Freeze a compact new-citation context review from two cached source pages."""
import datetime
import hashlib
import json
from pathlib import Path

WORK=Path('/workspace/research_work/mitl-upstream-20261019T0953')
OUT=WORK/'citation-review'
PDF=Path('/workspace/research_sources/ss304_properties/new_mills2004/mills_2004_official_pdf_response')
MANIFEST=Path('/workspace/Power-Flow/data/materials/mills_2004_source_manifest.json')

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def load(path):
    return json.loads(path.read_text())

status=load(Path('/workspace/research_work/mitl-upstream-status.json'))
assert status['active_lease_owner']=='mitl-upstream-20261019T0953'
assert datetime.datetime.now(datetime.timezone.utc)<datetime.datetime.fromisoformat(status['active_editing_lease_until_utc'])
assert sha(PDF)=='bf120dcbdf20d5f9b00494efd5e05643793edc0838deee8ba2c52988f9074836'
manifest=load(MANIFEST)
assert manifest['source']['doi']=='10.2355/isijinternational.44.1661'
assert manifest['source']['pdf_sha256']==sha(PDF)
inspection=load(WORK/'citation-inspect-two-pdf-pages.log')
assert inspection['pdf_sha256']==sha(PDF)
pages=inspection['pages']
assert [row['pdf_page_1based'] for row in pages]==[2,8]
context=pages[0]['matches'][0]['text']
reference_block=pages[1]['matches'][0]['text']
assert 'Density data for the liquid 316 alloy' in context
assert 'McCormick and Brooks8)' in context
assert '8)\nA. McCormick and R. F. Brooks:' in reference_block
assert 'PMP 1, 2, 3, Final Report, NPL, Teddington, Middlesex, UK,\n(1996).' in reference_block
record={
 'schema':'mitl-npl-processability-citation-context-review-v1',
 'status':'citation_identity_and_liquid316_context_confirmed',
 'reviewed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
 'scope':'Independent new-reference identity/material-use context from cached Mills2004 PDF pages2 and8 only; no repeated full source audit.',
 'source':{'title':manifest['source']['title'],'doi':manifest['source']['doi'],
   'cached_pdf':str(PDF),'pdf_sha256':sha(PDF),'source_manifest_sha256':sha(MANIFEST)},
 'new_citation':{'reference_number_in_mills':8,'authors':['A. McCormick','R. F. Brooks'],
   'literal_citation_normalized_whitespace':'A. McCormick and R. F. Brooks: MTS Programme on processability thermophysical property data for commercial alloys measured in PMP 1, 2, 3, Final Report, NPL, Teddington, Middlesex, UK, (1996).',
   'year':1996,'institution':'NPL, Teddington, Middlesex, UK',
   'report_number_or_DOI_supplied_by_this_citation':None,
   'identity_qualification':'The program/PMP wording and Final Report are citation facts; no unverified NPL report-series number, DOI or revised title is invented.'},
 'page_evidence':[{key:value for key,value in row.items() if key!='matches'} for row in pages],
 'pinpoints':[
  {'pdf_page_1based':2,'printed_page':1662,'location':'Left-column paragraph at the start of the density database discussion.',
   'source_block_bbox_pdf_points':pages[0]['matches'][0]['bbox'],
   'short_quote':'Density data for the liquid 316 alloy have been reported by McCormick and Brooks8)',
   'finding':'Mills specifically uses ref8 for liquid316 density. Its solid304/316 density lineage is refs6/7; liquid430 is ref9.'},
  {'pdf_page_1based':8,'printed_page':1668,'location':'REFERENCES, item8, right column, between Mills2002 item7 and Li et al. item9.',
   'containing_reference_block_bbox_pdf_points':pages[1]['matches'][0]['bbox'],
   'finding':'Bibliographic item8 identifies the1996 McCormick/Brooks NPL final report.'}],
 'material_use_qualification':{
  'liquid316_density_lineage_supported':True,'liquid304_density_established_by_citation':False,
  'all_grades_contained_in_new_report_known':False,'actual_dissertation_foil_properties_established':False,
  'interpretation':'The citation does not prove absence of304 from the underlying report. Actual grade/specimen/phase/composition/property coverage must be read from that report before any304 use.',
  'EOS_heat_capacity_or_phase_transition_input_adopted':False},
 'new_report_full_text_inspected_by_this_review':False,'network_or_source_requests':0,
 'scientific_model_or_numerical_tests_or_grids_executed':False,'physical_parameters_adopted':0,
 'public_files_changed':[],
 'source_containing_logs_excluded_from_delivery':[
  {'path':str(WORK/'citation-inspect-two-pdf-pages.log'),'sha256':sha(WORK/'citation-inspect-two-pdf-pages.log'),
   'reason':'Contains an extracted source reference block; retain locally and bind only its hash in review delivery.'},
  {'path':str(WORK/'citation-read-reference-context.log'),'sha256':sha(WORK/'citation-read-reference-context.log'),
   'reason':'Contains source text excerpts; retain locally and bind only its hash in review delivery.'}],
 'limitations':['Two-page citation/context verification only; source-author terminology retained.',
   'No liquid304 or actual-foil data inferred from a citation to liquid316.',
   'No same-alloy assumption, numerical parameter adoption or repeated Mills evaluation.'],
 'artifact_sha256':{'freeze_citation_context.py':sha(Path(__file__))}
}
path=OUT/'completion.json'
path.write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({'path':str(path),'sha256':sha(path),'status':record['status']}))
