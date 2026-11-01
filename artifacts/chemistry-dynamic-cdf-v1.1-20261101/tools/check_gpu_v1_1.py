import hashlib,json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE=Path('/home/treefor/projects/runtimes/pic-water-chemistry-v3.1.1-20261031')
ARR=Path('/home/treefor/projects/runtimes/pic-arrival-owner-v1-20261030/include/pf_arrival')
GPU=Path('/home/treefor/projects/runtimes/hall-gpu-resident-v1-20261030')
NVRTC=Path('/home/treefor/projects/runtimes/pic-cuda-gather-investigation-20261029/tools/nvrtc/nvidia/cuda_nvrtc')
CUDA=Path('/home/treefor/projects/runtimes/mfem-cuda-benchmark-20261029/cuda-libs/nvidia/cuda_runtime/include')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def verify(pins):
 for section in ('source_files','dependencies'):
  for p,h in pins[section].items():assert sha(ROOT/p)==h,p
 assert sha(Path(sys.executable).resolve())==pins['frozen_python_sha256']
def main():
 job=os.environ.get('SIMQ_JOB_ID');assert job and job.isdigit() and len(os.sched_getaffinity(0))==1 and os.environ.get('SIMQ_GPUS')
 assert not os.environ.get('PF_BUILD_ISA_FLAGS')
 pins=json.loads((ROOT/'GPU-V1.1-INPUTS.json').read_text());verify(pins);out=ROOT/'evidence'/('gpu-job-'+job);out.mkdir(exist_ok=False)
 def run(args,label):
  p=subprocess.run(args,cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=180);(out/(label+'.stdout')).write_bytes(p.stdout);(out/(label+'.stderr')).write_bytes(p.stderr);assert p.returncode==0,p.stderr.decode();return p.stdout
 run(['/usr/bin/c++','-std=c++20','-O2','-Wall','-Wextra','-Wno-c++26-extensions','-ffp-contract=off','-fno-fast-math','-I'+str(ROOT/'include/pf_chem_v3'),'-I'+str(ARR),'-I'+str(GPU),'-isystem',str(CUDA),'-isystem',str(NVRTC/'include'),str(ROOT/'tests/gpu_parity_v1_1.cpp'),str(GPU/'common/hall_gpu_common.cpp'),'-ldl','-lcrypto','-o',str(ROOT/'build/gpu_parity_v1_1')],'compile')
 run([str(ROOT/'build/gpu_parity_v1_1'),str(ROOT/'data/chemistry-cdf.blob'),str(out/'parity'),str(NVRTC/'lib/libnvrtc.so.12'),'GPU-a7b9b458-167f-cd57-fec9-3dea0995b6e9','2000','1',sha(ROOT/'data/chemistry-cdf.blob')],'parity')
 verify(pins);receipt={'schema':'chemistry-dynamic-cdf-gpu-receipt-v1','job':job,'allocated_cpu':sorted(os.sched_getaffinity(0)),'allocated_GPU':os.environ['SIMQ_GPUS'],'numerical_source_identity':pins['numerical_source_identity'],'inputs_sha256':sha(ROOT/'GPU-V1.1-INPUTS.json'),'binary_sha256':sha(ROOT/'build/gpu_parity_v1_1'),'files_sha256':{str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()},'passed':True,'production_qualified':False,'physical_model_validated':False,'native_adopted':False,'network_complete':False,'performance_benchmark':False};(out/'receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');print('PASS dynamic-CDF CPU/GPU populated events and unsupported selected-event refusal; shared GPU qualification, no timing claim')
if __name__=='__main__':main()
