#!/usr/bin/env python3
# burnside sum -> A048192 / A048193, checks the known values
import os
_HERE=os.path.dirname(os.path.abspath(__file__))
import importlib.util,sys,time
from math import factorial
from fractions import Fraction
def load(p,n):
    s=importlib.util.spec_from_file_location(n,os.path.join(_HERE,os.path.basename(p)));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
rd=load("rec_dp.py","rd")

def partitions(n,mx=None):
    if mx is None: mx=n
    if n==0: yield []; return
    for first in range(min(n,mx),0,-1):
        for rest in partitions(n-first,first):
            yield [first]+rest

def zlam(part):
    from collections import Counter
    c=Counter(part); z=1
    for i,m in c.items(): z*= (i**m)*factorial(m)
    return z

def counts_of(part):
    from collections import Counter
    return dict(Counter(part))

KNOWN_ALL ={1:1,2:2,3:4,4:10,5:27,6:94,7:393,8:2119,9:14524,10:126758,11:1392387,
            12:19109099,13:326005775,14:6905776799,15:181945055235}
KNOWN_CONN={1:1,2:1,3:2,4:5,5:15,6:58,7:272,8:1614,9:11911,10:109539,11:1247691,
            12:17566431,13:305310547,14:6558690953,15:174688164414}

def burnside(n):
    A=Fraction(0); C=Fraction(0)
    for part in partitions(n):
        cnt=counts_of(part)
        inst=rd.rd_inst(cnt) if hasattr(rd,'rd_inst') else rd.EquivDP.get(cnt.keys())
        z=zlam(part)
        fx=inst.fix_counts(cnt)
        cf=inst.cfix_tuple(tuple(cnt.get(s,0) for s in inst.S))
        A+=Fraction(fx,z); C+=Fraction(cf,z)
    assert A.denominator==1, f"A048193 n={n} not integer: {A}"
    assert C.denominator==1, f"A048192 n={n} not integer: {C}"
    return int(A),int(C)

if __name__=="__main__":
    NMAX=int(sys.argv[1]) if len(sys.argv)>1 else 16
    print(f"{'n':>3} {'A048193(all)':>16} {'A048192(conn)':>16}   check",flush=True)
    for n in range(1,NMAX+1):
        t0=time.time()
        a_all,a_conn=burnside(n)
        ck=""
        if n in KNOWN_ALL:
            ck=("all OK" if a_all==KNOWN_ALL[n] else f"all MISMATCH(exp {KNOWN_ALL[n]})")
            ck+=" | "+("conn OK" if a_conn==KNOWN_CONN[n] else f"conn MISMATCH(exp {KNOWN_CONN[n]})")
        print(f"{n:>3} {a_all:>16} {a_conn:>16}   {ck}  [{time.time()-t0:.1f}s]",flush=True)
