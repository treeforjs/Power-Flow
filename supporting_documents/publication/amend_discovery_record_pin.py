"""Add the missing identity pin for a preserved query-only discovery projection."""
from pathlib import Path
import hashlib,json
w=Path(__file__).resolve().parent
r=Path('/workspace/Power-Flow')
s=r/'supporting_documents'
def sha(data):return hashlib.sha256(data).hexdigest()
cat=s/'catalog.json'
old=cat.read_bytes()
(w/'catalog.before-discovery-record-pin.json').write_bytes(old)
d=json.loads(old)
p=s/'in_progress/mitl-upstream-20261019T0953/discovery_records.json'
b=p.read_bytes()
assert len(b)==1007 and sha(b)=='3f5c9bb45dd1fb54bfe16a609a628031eed0f673d1e451adfce93ce0b51339c6'
d['preserved_incomplete_work']['discovery_records']={'path':str(p.relative_to(r)),'bytes':len(b),'sha256':sha(b)}
cat.write_text(json.dumps(d,indent=2)+'\n')
local=w/'supporting-document-assembly.json'
public=s/'publication/assembly_record.json'
assert local.read_bytes()==public.read_bytes()
previous=local.read_bytes()
(w/'assembly.before-discovery-record-pin.json').write_bytes(previous)
a=json.loads(previous)
a['catalog']['sha256']=sha(cat.read_bytes())
a['catalog']['bytes']=cat.stat().st_size
updated=(json.dumps(a,indent=2)+'\n').encode()
local.write_bytes(updated)
public.write_bytes(updated)
record={'schema':'mitl-publication-discovery-record-pin-amendment-v1','status':'corrected_metadata_only','initial_verification_report_sha256':'c8c7f29169fe393638871021138bb71eceb9e212100aaad5bd9b35f4457a46e3','initial_issue':'The existing query-only discovery_records.json was omitted from the catalog identity list; independent subtree coverage rejected it.','catalog_before_sha256':sha(old),'catalog_after_sha256':sha(cat.read_bytes()),'assembly_before_sha256':sha(previous),'assembly_after_sha256':sha(updated),'added_pin':d['preserved_incomplete_work']['discovery_records'],'archive_records_unchanged':json.loads(old)['artifacts']==d['artifacts'],'existing_preserved_incomplete_fields_unchanged':all(d['preserved_incomplete_work'][k]==v for k,v in json.loads(old)['preserved_incomplete_work'].items()),'discovery_record_file_changed':False,'old_scientific_or_source_or_archive_producers_reexecuted':False,'resume_scope':'Failed and unreached publication metadata/link checks only; retain the initial successful archive readback and its frozen pins.'}
(s/'publication/discovery_record_pin_amendment.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
