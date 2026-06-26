#!/usr/bin/env python3
# all validation layers in one run
import os
_HERE=os.path.dirname(os.path.abspath(__file__))
import importlib.util,time
from collections import Counter
from math import factorial
from fractions import Fraction
def load(p,n):
    s=importlib.util.spec_from_file_location(n,os.path.join(_HERE,os.path.basename(p)));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
rd=load("rec_dp.py","rd")
bt=load("general_brute_bt.py","bt")
fails=0
def check(name,cond):
    global fails
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}",flush=True)
    if not cond: fails+=1

print("="*70)
print("LAYER 1: labeled limit  fix(1^n) == A058862 (all labeled chordal)")
A058862=[1,2,8,61,822,18154,617675,30888596,2192816760,215488096587]
lab=rd.EquivDP.get([1])
check("fix(1^n)=A058862, n=1..10", all(lab.fix((n,))==A058862[n-1] for n in range(1,11)))

print("="*70)
print("LAYER 2: blocks 1..10  cfix(orbit) == brute, every size incl composites")
b2=[{2:2,1:1},{3:2,1:1},{4:2,1:1},{5:1,2:1,1:1},{6:1,3:1},{6:2},{8:1,4:1},{8:1,2:1},
    {9:1,3:1},{10:1,5:1},{10:1,2:1},{4:1,2:1,1:1},{7:1,2:1},{4:2,2:1}]
ok=True
for cnt in b2:
    lam=[]; [lam.extend([s]*v) for s,v in cnt.items()]
    inst=rd.EquivDP.get(cnt.keys())
    dpcf=inst.cfix_tuple(tuple(cnt.get(s,0) for s in inst.S))
    _,cf,_=bt.solve(lam); ok&= (dpcf==cf)
check("cfix == brute for 14 mixed/composite cycle types (sizes 1..10)", ok)

print("="*70)
print("LAYER 3: fix(lambda) == brute (ALL invariant chordal), all blocks")
b3=[{2:1,1:1},{2:2},{1:3},{3:1,2:1},{4:1,1:1},{4:2},{6:1,2:1},{6:2},{8:1,2:1},
    {9:1,3:1},{10:1,2:1},{3:3},{2:4},{4:1,1:3},{5:1,1:1},{6:1,3:1}]
ok=True
for cnt in b3:
    lam=[]; [lam.extend([s]*v) for s,v in cnt.items()]
    inst=rd.EquivDP.get(cnt.keys()); dpfix=inst.fix_counts(cnt)
    fx,_,_=bt.solve(lam); ok&=(dpfix==fx)
check("fix == brute for 16 cycle types", ok)

print("="*70)
print("LAYER 4: every single orbit size 1..16  (fix,cfix) == brute")
ok=True
for c in range(1,17):
    inst=rd.EquivDP.get([c]); dpf=inst.fix((1,)); dpc=inst.cfix_tuple((1,))
    fx,cf,_=bt.solve([c]); ok&=(dpf==fx and dpc==cf)
check("single c-cycle fix,cfix == brute, c=1..16", ok)

print("="*70)
print("LAYER 5: Burnside reproduces all 15 KNOWN terms of A048192 & A048193")
def partitions(n,mx=None):
    if mx is None: mx=n
    if n==0: yield []; return
    for f in range(min(n,mx),0,-1):
        for r in partitions(n-f,f): yield [f]+r
def zlam(part):
    c=Counter(part); z=1
    for i,m in c.items(): z*=(i**m)*factorial(m)
    return z
def burnside(n):
    A=Fraction(0);C=Fraction(0)
    for part in partitions(n):
        cnt=dict(Counter(part)); inst=rd.EquivDP.get(cnt.keys()); z=zlam(part)
        A+=Fraction(inst.fix_counts(cnt),z)
        C+=Fraction(inst.cfix_tuple(tuple(cnt.get(s,0) for s in inst.S)),z)
    assert A.denominator==1 and C.denominator==1
    return int(A),int(C)
KALL=[1,2,4,10,27,94,393,2119,14524,126758,1392387,19109099,326005775,6905776799,181945055235]
KCON=[1,1,2,5,15,58,272,1614,11911,109539,1247691,17566431,305310547,6558690953,174688164414]
comp=[burnside(n) for n in range(1,16)]
check("A048193 (all) matches OEIS n=1..15",  all(comp[n-1][0]==KALL[n-1] for n in range(1,16)))
check("A048192 (conn) matches OEIS n=1..15", all(comp[n-1][1]==KCON[n-1] for n in range(1,16)))

print("="*70)
print("LAYER 6: NEW a(16) + Euler-transform self-consistency n=1..16")
a16=burnside(16)
print(f"    a(16): A048193(all)={a16[0]}  A048192(conn)={a16[1]}",flush=True)
check("a(16) all  == 5985406996403",  a16[0]==5985406996403)
check("a(16) conn == 5796153514484", a16[1]==5796153514484)
CONN=KCON+[a16[1]]; ALL=KALL+[a16[0]]
def euler(c):
    N=len(c); b=[0]*(N+1)
    for k in range(1,N+1):
        for m in range(k,N+1,k): b[m]+=k*c[k-1]
    a=[0]*(N+1); a[0]=1
    for n in range(1,N+1):
        s=b[n]+sum(b[k]*a[n-k] for k in range(1,n)); a[n]=s//n
    return a[1:]
check("Euler(A048192)==A048193 for all n=1..16 (incl a(16))", euler(CONN)==ALL)

print("="*70)
print(f"RESULT: {'ALL LAYERS PASS' if fails==0 else str(fails)+' FAILURES'}",flush=True)
