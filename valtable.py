#!/usr/bin/env python3
# per-cycle-type validation vs brute (prints latex rows)
import os
_HERE=os.path.dirname(os.path.abspath(__file__))
import importlib.util
from collections import Counter
def load(p,n):
    s=importlib.util.spec_from_file_location(n,os.path.join(_HERE,os.path.basename(p)));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
rd=load("rec_dp.py","rd")
bt=load("general_brute_bt.py","bt")
cts=[[1,1,1],[2,1,1],[2,2],[2,2,1],[3,1,1],[3,2,1],[3,3],[4,1,1],[4,2],[4,2,1],[4,4],
     [5,2],[5,2,1],[6,2],[6,3],[6,4],[8,2],[8,4],[9,3],[10,2],[10,5],[4,3,2],[3,3,3]]
allok=True; rows=[]
for ct in sorted(cts,key=lambda c:(sum(c),[-x for x in sorted(c)])):
    cnt=dict(Counter(ct)); n=sum(ct)
    inst=rd.EquivDP.get(cnt.keys())
    dpf=inst.fix_counts(cnt); dpc=inst.cfix_tuple(tuple(cnt.get(s,0) for s in inst.S))
    bf,bc,_=bt.solve(ct)
    ok=(dpf==bf and dpc==bc); allok&=ok
    parts=sorted(ct,reverse=True)
    # exponent notation like 4^2 2 1^3
    c=Counter(parts); lam=" ".join((f"{s}^{{{m}}}" if m>1 else f"{s}") for s,m in sorted(c.items(),reverse=True))
    tag=r"\checkmark" if ok else "FAIL"
    rows.append(f"${lam}$ & {n} & {dpf} & {dpc} & {tag} \\\\")
print("\n".join(rows))
print("% ALL MATCH" if allok else "% FAILURE")
