#!/usr/bin/env python3
# compute a(n) for a range and write the b-files
import os
_HERE=os.path.dirname(os.path.abspath(__file__))
import importlib.util,sys,time
from math import factorial
from fractions import Fraction
from collections import Counter
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
    c=Counter(part); z=1
    for i,m in c.items(): z*=(i**m)*factorial(m)
    return z
def burnside(n):
    A=Fraction(0);C=Fraction(0)
    for part in partitions(n):
        cnt=dict(Counter(part)); inst=rd.EquivDP.get(cnt.keys()); z=zlam(part)
        fx=inst.fix_counts(cnt); cf=inst.cfix_tuple(tuple(cnt.get(s,0) for s in inst.S))
        A+=Fraction(fx,z);C+=Fraction(cf,z)
    assert A.denominator==1 and C.denominator==1
    return int(A),int(C)

# known OEIS values, used to self-check n<=15
KALL ={1:1,2:2,3:4,4:10,5:27,6:94,7:393,8:2119,9:14524,10:126758,11:1392387,12:19109099,
       13:326005775,14:6905776799,15:181945055235}
KCON ={1:1,2:1,3:2,4:5,5:15,6:58,7:272,8:1614,9:11911,10:109539,11:1247691,12:17566431,
       13:305310547,14:6558690953,15:174688164414}

if __name__=="__main__":
    NMAX=int(sys.argv[1]) if len(sys.argv)>1 else 16
    allv={};connv={};bad=0
    fA=open("b048193_new.txt","w")
    fC=open("b048192_new.txt","w")
    for n in range(1,NMAX+1):
        t0=time.time(); a,c=burnside(n); allv[n]=a; connv[n]=c
        fA.write(f"{n} {a}\n"); fA.flush()
        fC.write(f"{n} {c}\n"); fC.flush()
        tag=""
        if n in KALL:
            ok=(a==KALL[n] and c==KCON[n]); bad+=0 if ok else 1
            tag="  known: OK" if ok else "  *** KNOWN MISMATCH ***"
        else:
            tag="  NEW"
        print(f"n={n:2d}  A048193(all)={a}  A048192(conn)={c}{tag}  [{time.time()-t0:.1f}s]",flush=True)
    fA.close(); fC.close()
    print("b-files written: b048193_new.txt, b048192_new.txt",flush=True)
    print("RESULT: all known terms (n<=15) reproduced" if bad==0
          else f"RESULT: {bad} known-term MISMATCHES",flush=True)
