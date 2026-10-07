"""Driver: rebuilds all regime certificates of Appendix B from the polynomial systems, without cached
solutions (branch_systems, p2_field, continuation of P1, P2 and P5, p1_slacks, endpoints)."""
import sys,subprocess,time,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
tasks=[['branch_systems.py'],['p2_field.py'],
       ['continue_certify.py','p1','4000'],
       ['continue_certify.py','p2','4000'],
       ['continue_certify.py','central','4000'],
       ['p1_slacks.py'],['endpoints.py']]
records=[]
for args in tasks:
    start=time.monotonic()
    print('START',*args,flush=True)
    proc=subprocess.run([sys.executable,'-u',*args],cwd=HERE)
    records.append(dict(command=args,exit_code=proc.returncode,seconds=round(time.monotonic()-start,3)))
    assert proc.returncode==0,args
print('ALL REGIME CERTIFICATES PASSED',flush=True)
print(json.dumps(records),flush=True)
