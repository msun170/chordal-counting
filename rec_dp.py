#!/usr/bin/env python3
# master equivariant dp: fix/cfix for any cycle type, via divisor-bundle recursion
import os
_HERE=os.path.dirname(os.path.abspath(__file__))
import importlib.util
from math import comb, gcd
from itertools import product
def load(p,n):
    s=importlib.util.spec_from_file_location(n,os.path.join(_HERE,os.path.basename(p)));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
asm=load("assemble.py","asm");asm.W=200

def divisors_gt1(c): return [d for d in range(2,c+1) if c%d==0]

class EquivDP:
    reg={}
    @classmethod
    def get(cls,sizes):
        S=tuple(sorted(set(sizes)))
        if S not in cls.reg: cls.reg[S]=EquivDP(S)
        return cls.reg[S]
    def __init__(self,S):
        self.S=S; self.n=len(S)
        self.idx={s:i for i,s in enumerate(S)}
        self._g1={};self._f={};self._g={};self._gT={};self._gH={};self._g2={};self._fT={};self._fH={}
    # tuple helpers
    def zero(self): return (0,)*self.n
    def add(self,a,b): return tuple(a[i]+b[i] for i in range(self.n))
    def sub(self,a,b): return tuple(a[i]-b[i] for i in range(self.n))
    def isz(self,a): return all(v==0 for v in a)
    def prodC(self,k,m):
        p=1
        for i in range(self.n): p*=comb(k[i],m[i])
        return p
    def subtuples(self,k):
        for m in product(*[range(k[i]+1) for i in range(self.n)]): yield m
    def markpos(self,k):
        for i in range(self.n-1,-1,-1):
            if k[i]>0: return i
        return -1
    def mark_coeff(self,k,kp,mp):
        cc=1
        for i in range(self.n):
            cc*= comb(k[i]-1,kp[i]-1) if i==mp else comb(k[i],kp[i])
        return cc
    # bundle GF machinery (multivariate over sub-types)
    def _subinfo(self,d):
        # returns (subS tuple, tau_of_size dict s->sub_type, g_of_size dict s->gcd(d,s))
        tau={s:s//gcd(d,s) for s in self.S}
        g={s:gcd(d,s) for s in self.S}
        subS=tuple(sorted(set(tau.values())))
        return subS,tau,g
    def _bundle_value(self,t,d,kp,x,z,kind,tshift=0):
        # peeled component-orbit is a d-fold bundle absorbing kp (all sizes with d|s). Returns
        # d^(K'-1) * sum_sel coeff(sel) * subInst.g1(t-tshift, sel-as-subboundary, k'-absorbed).
        subS,tau,g=self._subinfo(d)
        A=len(subS); sidx={s:i for i,s in enumerate(subS)}
        # absorbed sub-orbits k' (size s with d|s -> sub-type s//d == tau[s])
        kp_sub=[0]*A
        for i,s in enumerate(self.S):
            if kp[i]>0:
                if s%d!=0: return 0   # this bundle can't absorb size s
                kp_sub[sidx[tau[s]]]+=kp[i]
        kp_sub=tuple(kp_sub); Kp=sum(kp)
        if Kp==0: return 0
        # build coefficient poly over subS variables
        def emono(ti,j):
            key=[0]*A; key[sidx[ti]]=j; return tuple(key)
        def mmul(a,b):
            r={}
            for ka,va in a.items():
                for kb,vb in b.items():
                    kk=tuple(ka[i]+kb[i] for i in range(A))
                    r[kk]=r.get(kk,0)+va*vb
            return r
        def full_from(counts):  # (1+v_tau)^{g_s*counts[s]} product
            poly={(0,)*A:1}
            for s in self.S:
                p=g[s]*counts[self.idx[s]]
                if p:
                    fac={emono(tau[s],j):comb(p,j) for j in range(p+1)}
                    poly=mmul(poly,fac)
            return poly
        def cover_from(counts):  # per-orbit nonempty: ((1+v)^{g_s}-1)^{counts[s]}
            poly={(0,)*A:1}
            for s in self.S:
                cnt=counts[self.idx[s]]
                if cnt:
                    base={emono(tau[s],j):comb(g[s],j) for j in range(1,g[s]+1)}  # j>=1
                    fac={(0,)*A:1}
                    for _ in range(cnt): fac=mmul(fac,base)
                    poly=mmul(poly,fac)
            return poly
        full=full_from(x)
        if kind=='all':
            fz=full_from(z); diff={k:full.get(k,0)-fz.get(k,0) for k in set(full)|set(fz)}
        elif kind=='cover':
            diff=cover_from(x)
        else: # notall
            if z==x: return 0
            fz=full_from(z); cov=cover_from(x)
            keys=set(full)|set(fz)|set(cov)
            diff={k:full.get(k,0)-fz.get(k,0)-cov.get(k,0) for k in keys}
        sub=EquivDP.get(subS)
        s=0
        for sel,co in diff.items():
            if co==0: continue
            if all(v==0 for v in sel) and kind!='cover': continue
            s+=co*sub.g1(t-tshift,sel,kp_sub)
        return (d**(Kp-1))*s
    # core recurrences (mirror validated orbit_dp/block124, multi-divisor bundles)
    def g1(self,t,x,k):
        if t==0 or self.isz(k): return 0
        key=(t,x,k)
        if key in self._g1: return self._g1[key]
        s=0
        for l in self.subtuples(k):
            if self.isz(l): continue
            s+=self.prodC(k,l)*self.f(t,x,l,self.sub(k,l))
        self._g1[key]=s;return s
    def f(self,t,x,l,k):
        if t==0: return 0
        if t==1: return 1 if self.isz(k) else 0
        if self.isz(k): return 0
        key=(t,x,l,k)
        if key in self._f: return self._f[key]
        s=0; xl=self.add(x,l)
        for m in self.subtuples(k):
            s+=self.prodC(k,m)*self.fT(t,x,l,m)*self.g(t-2,xl,self.sub(k,m),x)
        self._f[key]=s;return s
    def g(self,t,x,k,z):
        if t==0: return 1 if self.isz(k) else 0
        if t<0: return 0
        key=(t,x,k,z)
        if key in self._g: return self._g[key]
        s=0
        for m in self.subtuples(k):
            s+=self.prodC(k,m)*self.gT(t,x,m,z)*self.g(t-1,x,self.sub(k,m),z)
        self._g[key]=s;return s
    def _Atau(self,t,x,z,kp,notall):
        s=0
        for xp in self.subtuples(x):
            if self.isz(xp): continue
            if notall and xp==x: continue
            co=self.prodC(x,xp)-self.prodC(z,xp)
            if co: s+=co*self.g1(t,xp,kp)
        return s
    def gT(self,t,x,k,z):
        if self.isz(k): return 1
        key=(t,x,k,z)
        if key in self._gT: return self._gT[key]
        mp=self.markpos(k); c=self.S[mp]; tot=0
        for kp in self.subtuples(k):
            if kp[mp]<1: continue
            cc=self.mark_coeff(k,kp,mp); rest=self.sub(k,kp)
            term=self._Atau(t,x,z,kp,False)
            for d in divisors_gt1(c):
                term+=self._bundle_value(t,d,kp,x,z,'all')
            tot+=cc*term*self.gT(t,x,rest,z)
        self._gT[key]=tot;return tot
    def gH(self,t,x,k,z):
        if self.isz(k): return 1
        if t<1: return 0
        key=(t,x,k,z)
        if key in self._gH: return self._gH[key]
        mp=self.markpos(k); c=self.S[mp]; tot=0
        for kp in self.subtuples(k):
            if kp[mp]<1: continue
            cc=self.mark_coeff(k,kp,mp); rest=self.sub(k,kp)
            term=self._Atau(t,x,z,kp,True)
            for d in divisors_gt1(c):
                term+=self._bundle_value(t,d,kp,x,z,'notall')
            tot+=cc*term*self.gH(t,x,rest,z)
        self._gH[key]=tot;return tot
    def g2(self,t,x,k):
        if self.isz(k) or t==0: return 0
        key=(t,x,k)
        if key in self._g2: return self._g2[key]
        def R(m): return self.g1(t,x,m)+self.g2(t,x,m)
        def Rest0(m): return (1 if self.isz(m) else 0)+R(m)
        mp=self.markpos(k); c=self.S[mp]; tot=0
        for kp in self.subtuples(k):
            if kp[mp]<1: continue
            cc=self.mark_coeff(k,kp,mp); rest=self.sub(k,kp)
            tot+=cc*self.g1(t,x,kp)*R(rest)
            bcov=0
            for d in divisors_gt1(c):
                bcov+=self._bundle_value(t,d,kp,x,self.zero(),'cover')
            if bcov: tot+=cc*bcov*Rest0(rest)
        self._g2[key]=tot;return tot
    def fT(self,t,x,l,k):
        if t==0 or t==1 or self.isz(k): return 0
        key=(t,x,l,k)
        if key in self._fT: return self._fT[key]
        xl=self.add(x,l)
        s=self.fH(t,x,x,l,k)
        for m in self.subtuples(k):
            if self.isz(m): continue
            s+=self.prodC(k,m)*(self.g1(t-1,xl,m)*self.fH(t,x,x,l,self.sub(k,m))
                               +self.g2(t-1,xl,m)*self.gH(t-1,xl,self.sub(k,m),x))
        self._fT[key]=s;return s
    def fH(self,t,x,z,l,k):
        if t<2 or self.isz(k): return 0
        key=(t,x,z,l,k)
        if key in self._fH: return self._fH[key]
        mp=self.markpos(k); c=self.S[mp]; tot=0
        def weight(xp,lp):
            anyL=not self.isz(lp)
            base=self.prodC(x,xp)
            if not anyL: base-=self.prodC(z,xp)
            return base*self.prodC(l,lp)
        def nxt(lp,kp,rest):
            if lp==l: return self.gH(t-1,self.add(x,lp),rest,z)
            return self.fH(t,self.add(x,lp),z,self.sub(l,lp),rest)
        for kp in self.subtuples(k):
            if kp[mp]<1: continue
            cc=self.mark_coeff(k,kp,mp); rest=self.sub(k,kp)
            # d=1 TypeA: recursive g1 on touched orbits
            sA=0
            for xp in self.subtuples(x):
                for lp in self.subtuples(l):
                    if self.isz(xp) and self.isz(lp): continue
                    if xp==x and lp==l: continue
                    w=weight(xp,lp)
                    if w==0: continue
                    gg=self.g1(t-1,self.add(xp,lp),kp)
                    if gg==0: continue
                    sA+=w*gg*nxt(lp,kp,rest)
            tot+=cc*sA
            # bundles d>1 (TypeB): cover_touched per touched orbit, sub-instance g1
            for d in divisors_gt1(c):
                if any(kp[i]>0 and self.S[i]%d!=0 for i in range(self.n)): continue
                sB=0
                for xp in self.subtuples(x):
                    for lp in self.subtuples(l):
                        if self.isz(xp) and self.isz(lp): continue
                        if xp==x and lp==l: continue
                        w=weight(xp,lp)
                        if w==0: continue
                        seen=self.add(xp,lp)
                        gv=self._bundle_touched(t,d,kp,seen)
                        if gv==0: continue
                        sB+=w*gv*nxt(lp,kp,rest)
                tot+=cc*sB
        self._fH[key]=tot;return tot
    def _bundle_touched(self,t,d,kp,seen):
        # absorption bundle: touched orbits each nonempty (cover_touched), sub-instance g1 at t-1
        subS,tau,g=self._subinfo(d)
        A=len(subS); sidx={s:i for i,s in enumerate(subS)}
        kp_sub=[0]*A
        for i,s in enumerate(self.S):
            if kp[i]>0:
                if s%d!=0: return 0
                kp_sub[sidx[tau[s]]]+=kp[i]
        kp_sub=tuple(kp_sub); Kp=sum(kp)
        if Kp==0: return 0
        def emono(ti,j):
            key=[0]*A; key[sidx[ti]]=j; return tuple(key)
        def mmul(a,b):
            r={}
            for ka,va in a.items():
                for kb,vb in b.items():
                    kk=tuple(ka[i]+kb[i] for i in range(A))
                    r[kk]=r.get(kk,0)+va*vb
            return r
        poly={(0,)*A:1}
        for i,s in enumerate(self.S):
            cnt=seen[i]
            if cnt:
                base={emono(tau[s],j):comb(g[s],j) for j in range(1,g[s]+1)}
                fac={(0,)*A:1}
                for _ in range(cnt): fac=mmul(fac,base)
                poly=mmul(poly,fac)
        sub=EquivDP.get(subS)
        s=0
        for sel,co in poly.items():
            if co==0 or all(v==0 for v in sel): continue
            s+=co*sub.g1(t-1,sel,kp_sub)
        return (d**(Kp-1))*s
    def c_orbit(self,counts):
        k=tuple(counts.get(s,0) for s in self.S)
        n=sum(s*counts.get(s,0) for s in self.S)
        return sum(self.g1(t,self.zero(),k) for t in range(1,n+2))
    # ALL invariant chordal (fix) via component-bundle decomposition
    def nvert(self,k): return sum(self.S[i]*k[i] for i in range(self.n))
    def cfix_tuple(self,k):
        if self.isz(k): return 1
        return sum(self.g1(t,self.zero(),k) for t in range(1,self.nvert(k)+2))
    def bundle(self,mu):
        # # single pi-orbit-of-components on orbit-multiset mu (d=1 connected + d>1 isomorphic bundles)
        tot=self.cfix_tuple(mu)
        present=[self.S[i] for i in range(self.n) if mu[i]>0]
        if not present: return tot
        G=0
        for s in present: G=gcd(G,s)
        M=sum(mu)
        for d in divisors_gt1(G):
            sub_counts={}
            for i in range(self.n):
                if mu[i]>0:
                    st=self.S[i]//d; sub_counts[st]=sub_counts.get(st,0)+mu[i]
            subInst=EquivDP.get(sub_counts.keys())
            smu=tuple(sub_counts.get(s,0) for s in subInst.S)
            tot+=(d**(M-1))*subInst.cfix_tuple(smu)
        return tot
    def fix(self,k):
        if self.isz(k): return 1
        if not hasattr(self,'_fix'): self._fix={}
        if k in self._fix: return self._fix[k]
        mp=self.markpos(k); tot=0
        for mu in self.subtuples(k):
            if mu[mp]<1: continue
            cc=self.mark_coeff(k,mu,mp)
            tot+=cc*self.bundle(mu)*self.fix(self.sub(k,mu))
        self._fix[k]=tot;return tot
    def fix_counts(self,counts):
        return self.fix(tuple(counts.get(s,0) for s in self.S))

def c_orbit(counts):
    return EquivDP.get(counts.keys()).c_orbit(counts)

if __name__=="__main__":
    import sys
    # Stage 1: S={1} reproduces labeled A007134
    lab=EquivDP.get([1]); A007134=[1,1,4,35,541,13302,489287]
    got=[sum(lab.g1(t,(0,),(n,)) for t in range(1,n+2)) for n in range(1,8)]
    print("S={1} labeled:",got,"->","MATCH" if got==A007134 else "MISMATCH",flush=True)
    # Stage 2: S={1,2} reproduces cfix(2^a 1^b)
    d2=EquivDP.get([1,2])
    exp2={(1,1):2,(1,2):7,(2,1):25,(2,2):158,(3,1):817}
    ok=all(d2.c_orbit({2:a,1:b})==e for (a,b),e in exp2.items())
    print("S={1,2} involution cfix:","MATCH" if ok else "MISMATCH",flush=True)
