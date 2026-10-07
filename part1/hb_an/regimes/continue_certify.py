"""Certified continuation over 4000 adjacent parameter slabs (Table 2).

Usage: python continue_certify.py p1|p2|central [n]. For every slab the Krawczyk test of interval_core
is applied, the slab checks of branch_checks are evaluated, and an enclosure obtained at the right end of
a slab is required to lie in the box of the next slab, so that the certified zeros form one branch.
Output: <name>_boxes.json."""
import sys,pickle,json,numpy as np
from pathlib import Path
from interval_core import run
from branch_checks import check_p1,check_p2,check_central
name=sys.argv[1];n=int(sys.argv[2]) if len(sys.argv)>2 else 1000
par,vs,F,Q,su=pickle.loads(Path(name+'_system.pkl').read_bytes())
if name=='p1':grid=np.linspace(0,.733,n+1);guess=[.125,-.25]
elif name=='p2':grid=np.linspace(1.02358,1.03330,n+1);guess=[.51156,-.5038,.2797,.0516,.04556,.36489]
else:grid=np.linspace(0,.023582,n+1);guess=[.51466,-.25854,.23189,2.579,2.1051,.374348,.596502]
check={"p1":lambda:check_p1(Q),"p2":lambda:check_p2(vs),"central":lambda:check_central(Q)}[name]()
out=run(F,par,vs,grid,guess,check=check,name=name)
Path(name+'_boxes.json').write_text(json.dumps(out,indent=2))
print(name,'PASS',len(out),'boxes')
