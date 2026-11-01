import importlib.util, pathlib, tempfile
ROOT=pathlib.Path(__file__).resolve().parents[1]
s=importlib.util.spec_from_file_location('primary_cdf',ROOT/'src/primary_cdf.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
raw=pathlib.Path('/home/treefor/projects/runtimes/water-fragment-source-inventory-20261101/references/sartori2024-table-I.txt').read_bytes()
assert len(m.parse(raw))==100
bad=[b'',raw[:20],raw+b'x'*65537,raw.replace(b'b3S_u+',b'a3S_g+'),raw.replace(b'v=0',b'v=1'),raw.replace(b'0.0010000',b'0.0000000'),raw.replace(b'0.999\t',b'1.000\t'),raw.replace(b'1.342054',b'NaN'),raw.replace(b'1.342054',b'-1.342054'),raw.replace(b'1.342054',b'1.342054\t1'),raw.replace(b'1.639446',b'1.342054'),raw.replace(b'0.01108081',b'0.0010000'),raw+b'\xff']
for value in bad:
    try:m.parse(value)
    except (ValueError,UnicodeError):pass
    else:raise AssertionError('malformed source accepted')
with tempfile.TemporaryDirectory() as tmp:
    p=pathlib.Path(tmp)/'source';p.write_bytes(raw+b'\n')
    try:m.load(p)
    except ValueError:pass
    else:raise AssertionError('changed identity accepted')
print('PASS 14 malformed/identity controls; exact 100 primary knots')
