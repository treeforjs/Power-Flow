"""First OCR of the acquired scanned NPL density-method report; raw pages stay cached."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import datetime,hashlib,json,subprocess,time
import fitz
P=Path('/workspace/research_sources/npl_processability_1996/density_method_pdf.body');C=Path('/workspace/research_sources/npl_processability_1996/density-method-pages');W=Path(__file__).resolve().parent
expected='833b5680249c93b0f903e5143f9e5bbef62677d92f11ad416bda6866bdce5fda'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(P)==expected
version=subprocess.run(['tesseract','--version'],check=True,capture_output=True,text=True)
(C/'tesseract-version.txt').write_text(version.stdout+version.stderr)
doc=fitz.open(P);assert len(doc)==25
pages=[]
for i,page in enumerate(doc):
 image=C/f'page-{i+1:02d}.png'
 if image.exists():raise RuntimeError('fresh OCR render would overwrite cache')
 pix=page.get_pixmap(dpi=250);pix.save(image)
 pages.append({'pdf_page_1based':i+1,'pdf_page_index':i,'pdf_rect_pt':list(page.rect),'render_dpi':250,'image':str(image),'image_pixels':[pix.width,pix.height],'image_sha256':sha(image)})
def ocr(r):
 image=Path(r['image']);out=C/image.stem
 command=['tesseract',str(image),str(out),'--dpi','250','--psm','3','txt','tsv']
 start=datetime.datetime.now(datetime.timezone.utc).isoformat();t=time.monotonic();p=subprocess.run(command,capture_output=True,text=True)
 (C/(image.stem+'.stderr')).write_text(p.stderr)
 return {**r,'OCR_command':command,'OCR_start_utc':start,'OCR_elapsed_s':time.monotonic()-t,'OCR_exit_code':p.returncode,'OCR_text_path':str(out.with_suffix('.txt')),'OCR_text_sha256':sha(out.with_suffix('.txt')) if out.with_suffix('.txt').exists() else None,'OCR_tsv_sha256':sha(out.with_suffix('.tsv')) if out.with_suffix('.tsv').exists() else None,'OCR_stderr_sha256':sha(C/(image.stem+'.stderr'))}
with ThreadPoolExecutor(max_workers=4) as ex:results=list(ex.map(ocr,pages))
assert all(r['OCR_exit_code']==0 for r in results)
report={'schema':'npl-density-method-first-OCR-v1','status':'all_25_pages_rendered_and_OCR_completed','source_pdf_path':str(P),'source_pdf_sha256':expected,'source_pdf_bytes':P.stat().st_size,'pdf_pages':len(doc),'render_library':'PyMuPDF '+fitz.VersionBind,'tesseract_version_file_sha256':sha(C/'tesseract-version.txt'),'raw_cache_directory':str(C),'raw_text_images_are_not_for_archive':True,'pages':results,'script_sha256':sha(Path(__file__)),'scientific_property_evaluation_or_fit':False,'external_requests':0}
(W/'ocr_inventory.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'status':report['status'],'pages':len(results),'source_pdf_sha256':expected,'inventory_sha256':sha(W/'ocr_inventory.json')},indent=2))
