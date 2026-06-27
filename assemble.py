#!/usr/bin/env python3
# labeled hlv chordal dp (the base case) + the involution tau dp
import os
_HERE=os.path.dirname(os.path.abspath(__file__))
import sys
from math import comb
from functools import lru_cache

W = 0  # clique bound (pairs for tau / vertices for labeled), set per target

# labeled HLV DP
import functools
def memo(f):
    c={}
    @functools.wraps(f)
    def w(*a):
        if a in c: return c[a]
        r=f(*a); c[a]=r; return r
    return w

@memo
def g_l(t,x,z,k):
    if t==0: return 1 if k==0 else 0
    s=0
    for kk in range(0,k+1):
        s+=comb(k,kk)*gT_l(t,x,z,kk)*g_l(t-1,x,z,k-kk)
    return s
@memo
def gT_l(t,x,z,k):
    if k==0: return 1
    s=0
    for kk in range(1,k+1):
        for xx in range(1,x+1):
            s+=comb(k-1,kk-1)*(comb(x,xx)-comb(z,xx))*g1_l(t,xx,kk)*gT_l(t,x,z,k-kk)
    return s
@memo
def gH_l(t,x,z,k):
    if k==0: return 1
    s=0
    for kk in range(1,k+1):
        for xx in range(1,x):
            s+=comb(k-1,kk-1)*(comb(x,xx)-comb(z,xx))*g1_l(t,xx,kk)*gH_l(t,x,z,k-kk)
    return s
@memo
def g1_l(t,x,k):
    if k==0 or t==0: return 0
    s=0
    for l in range(1,k+1):
        s+=comb(k,l)*f_l(t,x,l,k-l)
    return s
@memo
def g2_l(t,x,k):
    if k==0 or t==0: return 0
    s=0
    for kk in range(1,k):
        s+=comb(k-1,kk-1)*g1_l(t,x,kk)*(g1_l(t,x,k-kk)+g2_l(t,x,k-kk))
    return s
@memo
def f_l(t,x,l,k):
    if x+l>W: return 0
    if t==0: return 0
    if t==1: return 1 if k==0 else 0
    if k==0: return 0
    s=0
    for kk in range(1,k+1):
        s+=comb(k,kk)*fT_l(t,x,l,kk)*g_l(t-2,x+l,x,k-kk)
    return s
@memo
def fT_l(t,x,l,k):
    if t==0 or t==1 or k==0: return 0
    s=fH_l(t,x,x,l,k)
    for kk in range(1,k):
        s+=comb(k,kk)*g1_l(t-1,x+l,kk)*fH_l(t,x,x,l,k-kk)
    for kk in range(1,k+1):
        s+=comb(k,kk)*g2_l(t-1,x+l,kk)*gH_l(t-1,x+l,x,k-kk)
    return s
@memo
def fH_l(t,x,z,l,k):
    if t==0 or t==1 or k==0: return 0
    s=0
    for kk in range(1,k+1):
        for xx in range(0,x+1):
            for ll in range(0,l+1):
                if xx+ll==0 or xx+ll==x+l: continue
                co=comb(x,xx) if ll>0 else comb(x,xx)-comb(z,xx)
                nxt = fH_l(t,x+ll,z,l-ll,k-kk) if ll<l else gH_l(t-1,x+ll,z,k-kk)
                s+=comb(k-1,kk-1)*comb(l,ll)*co*g1_l(t-1,xx+ll,kk)*nxt
    return s
def chordal_conn_l(k):
    return sum(g1_l(t,0,k) for t in range(1,k+1))

# tau equivariant DP, all 8 validated recurrences
# coefficients
@memo
def Mc(x,xv):  # #xv-subsets of 2x verts covering all x pairs
    return sum((-1)**i*comb(x,i)*comb(2*x-2*i,xv) for i in range(0,x+1))
@memo
def Mp(x,z,xv):  # #xv-subsets missing >=1 pair AND beyond first 2z verts
    miss=comb(2*x,xv)-Mc(x,xv)
    if z>=x: return 0
    return miss-comb(2*z,xv)
@memo
def Gp(p,kk,t):  # labeled comp seeing p boundary PAIRS via a/b/both modes
    return sum(comb(p,v-p)*(2**(2*p-v))*g1_l(t,v,kk) for v in range(p,2*p+1))

@memo
def g1_t(t,x,k):
    if k==0 or t==0: return 0
    return sum(comb(k,l)*f_t(t,x,l,k-l) for l in range(1,k+1))
@memo
def f_t(t,x,l,k):
    if x+l>W: return 0
    if t==0: return 0
    if t==1: return 1 if k==0 else 0
    if k==0: return 0
    return sum(comb(k,kk)*fT_t(t,x,l,kk)*g_t(t-2,x+l,x,k-kk) for kk in range(1,k+1))
@memo
def g_t(t,x,z,k):
    if t==0: return 1 if k==0 else 0
    return sum(comb(k,kk)*gT_t(t,x,z,kk)*g_t(t-1,x,z,k-kk) for kk in range(0,k+1))
@memo
def gT_t(t,x,z,k):
    if k==0: return 1
    tot=0
    for kk in range(1,k+1):
        A=sum((comb(x,xx)-comb(z,xx))*g1_t(t,xx,kk) for xx in range(1,x+1))
        B=sum((comb(2*x,xv)-comb(2*z,xv))*(2**(kk-1))*g1_l(t,xv,kk) for xv in range(1,2*x+1))
        tot+=comb(k-1,kk-1)*(A+B)*gT_t(t,x,z,k-kk)
    return tot
@memo
def gH_t(t,x,z,k):
    if k==0: return 1
    tot=0
    for kk in range(1,k+1):
        A=sum((comb(x,xx)-comb(z,xx))*g1_t(t,xx,kk) for xx in range(1,x))
        B=sum(Mp(x,z,xv)*(2**(kk-1))*g1_l(t,xv,kk) for xv in range(1,2*x+1))
        tot+=comb(k-1,kk-1)*(A+B)*gH_t(t,x,z,k-kk)
    return tot
@memo
def g2_t(t,x,k):
    if k==0 or t==0: return 0
    def R(m): return g1_t(t,x,m)+g2_t(t,x,m)
    def Rest0(m): return (1 if m==0 else 0)+R(m)
    tot=0
    for kk in range(1,k+1):
        Bpeel=(2**(kk-1))*sum(Mc(x,xv)*g1_l(t,xv,kk) for xv in range(x,2*x+1)) if x>0 else (2**(kk-1))*g1_l(t,0,kk)
        tot+=comb(k-1,kk-1)*( g1_t(t,x,kk)*R(k-kk) + Bpeel*Rest0(k-kk) )
    return tot
@memo
def fT_t(t,x,l,k):
    if t==0 or t==1 or k==0: return 0
    s=fH_t(t,x,x,l,k)
    for kk in range(1,k):
        s+=comb(k,kk)*g1_t(t-1,x+l,kk)*fH_t(t,x,x,l,k-kk)
    for kk in range(1,k+1):
        s+=comb(k,kk)*g2_t(t-1,x+l,kk)*gH_t(t-1,x+l,x,k-kk)
    return s
@memo
def fH_t(t,x,z,l,k):
    if t==0 or t==1 or k==0: return 0
    tot=0
    for kk in range(1,k+1):
        cb=comb(k-1,kk-1)
        # TypeA: tau-inv component sees xx X-pairs + ll L-pairs
        for xx in range(0,x+1):
            for ll in range(0,l+1):
                if xx+ll==0 or xx+ll==x+l: continue
                co=comb(x,xx) if ll>0 else comb(x,xx)-comb(z,xx)
                nxt = fH_t(t,x+ll,z,l-ll,k-kk) if ll<l else gH_t(t-1,x+ll,z,k-kk)
                tot+=cb*comb(l,ll)*co*g1_t(t-1,xx+ll,kk)*nxt
        # TypeB: swapped pair, half sees xp X-pairs + lp L-pairs (modes via Gp)
        for xp in range(0,x+1):
            for lp in range(0,l+1):
                if xp+lp==0 or xp+lp==x+l: continue
                co=comb(x,xp) if lp>0 else comb(x,xp)-comb(z,xp)
                nxt = fH_t(t,x+lp,z,l-lp,k-kk) if lp<l else gH_t(t-1,x+lp,z,k-kk)
                tot+=cb*(2**(kk-1))*co*comb(l,lp)*Gp(xp+lp,kk,t-1)*nxt
    return tot
def cfix_t(j):
    return sum(g1_t(t,0,j) for t in range(1,j+1))

if __name__=="__main__":
    mode=sys.argv[1] if len(sys.argv)>1 else "both"
    if mode in ("lab","both"):
        W=8
        A007134=[1,1,4,35,541,13302,489287,25864897]
        got=[chordal_conn_l(n) for n in range(1,9)]
        print("labeled:", got, "->", "MATCH" if got==A007134 else "MISMATCH")
    if mode in ("tau","both"):
        W=100
        cfix=[1,7,154,7289,645981,100906896]
        got=[cfix_t(j) for j in range(1,7)]
        print("cfix  got:   ", got)
        print("cfix expect: ", cfix)
        print("==> MATCH!!" if got==cfix else "==> MISMATCH")
    if mode=="compute":
        import time
        W=100
        N=int(sys.argv[2]) if len(sys.argv)>2 else 9
        # cfix(j) = connected tau-inv chordal on j pairs; c_j = labeled connected (A007134)
        cf=[]; cj=[]
        for j in range(1,N+1):
            t0=time.time(); v=cfix_t(j); cf.append(v); cj.append(chordal_conn_l(j))
            print(f"cfix({j}) = {v}    [{time.time()-t0:.1f}s]", flush=True)
        # fix(m) = ALL tau-inv chordal on m pairs via connected reduction (validated)
        fix=[1]  # fix(0)=1
        for m in range(1,N+1):
            s=0
            for j in range(1,m+1):
                s+=comb(m-1,j-1)*(cf[j-1]+(2**(j-1))*cj[j-1])*fix[m-j]
            fix.append(s)
        print("fix(2^m) m=1..:", fix[1:])
