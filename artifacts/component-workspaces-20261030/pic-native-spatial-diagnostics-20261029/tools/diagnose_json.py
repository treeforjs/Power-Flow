#!/usr/bin/env python3
import argparse,pathlib,sys,json,hashlib
P=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'source'))
from diagnostics import diagnose
def main():
 p=argparse.ArgumentParser();p.add_argument('--input',type=pathlib.Path,required=True);p.add_argument('--out',type=pathlib.Path,required=True);a=p.parse_args();assert not a.out.exists();assert a.input.stat().st_size<=128*1024**2
 raw=a.input.read_bytes();data=json.loads(raw);result=diagnose(data);result['input_file_sha256']=hashlib.sha256(raw).hexdigest();text=json.dumps(result,indent=2,allow_nan=False)+'\n';assert len(text.encode())<=128*1024**2;a.out.write_text(text)
 print('Read-only spatial/energy diagnostics complete; no evolution or full-mechanism qualification')
if __name__=='__main__':main()
