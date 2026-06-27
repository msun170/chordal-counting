#!/usr/bin/env python3
# machine checks for the structural lemmas
import os
_HERE=os.path.dirname(os.path.abspath(__file__))
import importlib.util
from math import gcd
from itertools import combinations
def load(p,n):
    s=importlib.util.spec_from_file_location(n,os.path.join(_HERE,os.path.basename(p)));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
bt=load("general_brute_bt.py","bt")

def lcm(a,b): return a*b//gcd(a,b)
def perm_of(ct):
    cyc,n=bt.cycles_of(ct); pi=[0]*n
    for C in cyc:
        for i in range(len(C)): pi[C[i]]=C[(i+1)%len(C)]
    return pi,cyc,n
def all_invariant_graphs(ct):
    cyc,n=bt.cycles_of(ct)
    orbs=[]
    for k in range(len(cyc)):
        orbs+=bt.within_orbits(cyc[k])
        for j in range(k): orbs+=bt.between_orbits(cyc[k],cyc[j])
    for mask in range(1<<len(orbs)):
        adj=[set() for _ in range(n)]
        mm=mask;b=0
        while mm:
            if mm&1:
                for (u,v) in orbs[b]: adj[u].add(v);adj[v].add(u)
            mm>>=1;b+=1
        yield adj,n
def is_chordal(adj,n): return bt.is_chordal_set(adj,set(range(n)))
def evap_times(adj,n):
    alive=set(range(n)); t=0; tm={}
    while alive:
        t+=1; simp=[]
        for v in alive:
            Nv=(adj[v]&alive)-{v}; ok=True
            for u in Nv:
                if (adj[u]&Nv)!=(Nv-{u}): ok=False;break
            if ok: simp.append(v)
        if not simp: return None   # not chordal
        for v in simp: tm[v]=t
        alive-=set(simp)
    return tm
def components(adj,n):
    seen=set(); comps=[]
    for s in range(n):
        if s in seen: continue
        st=[s];seen.add(s);cm={s}
        while st:
            v=st.pop()
            for u in adj[v]:
                if u not in seen: seen.add(u);st.append(u);cm.add(u)
        comps.append(frozenset(cm))
    return comps
def apply_pi_set(pi,S): return frozenset(pi[v] for v in S)

fails=0
def report(name,ok):
    global fails
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}",flush=True)
    if not ok: fails+=1

#  L3: pure group theory
ok=True
for s in range(1,13):
    for d in range(1,13):
        # pi = +1 on Z/s ; pi^d = +d ; orbits of <d> in Z/s
        seen=set();norb=0;sizes=set()
        for start in range(s):
            if start in seen: continue
            cur=start;sz=0
            while cur not in seen:
                seen.add(cur);cur=(cur+d)%s;sz+=1
            norb+=1;sizes.add(sz)
        if norb!=gcd(d,s): ok=False
        if sizes and (sizes!={s//gcd(d,s)}): ok=False
report("L3 block decomposition: #orbits=gcd(d,s), size=s/gcd(d,s)  (all 1<=d,s<=12)",ok)

# L1 & L2: over enumerated invariant chordal graphs
types=[[2,1],[2,2],[3,1],[4,1],[3,2],[2,2,1],[4,2],[6],[2,1,1],[3,3],[4,1,1]]
L1ok=True;L2ok=True;count=0
for ct in types:
    pi,cyc,n=perm_of(ct); N=1
    for c in ct: N=lcm(N,c)
    for adj,_ in all_invariant_graphs(ct):
        tm=evap_times(adj,n)
        if tm is None: continue   # non-chordal
        if not is_chordal(adj,n): continue
        count+=1
        # L1: orbit-mates share evap time
        for C in cyc:
            if len({tm[v] for v in C})!=1: L1ok=False
        # L2: component orbit size divides N, component is <pi^d>-invariant
        comps=set(components(adj,n))
        for K in comps:
            orb=[K];cur=apply_pi_set(pi,K)
            while cur!=K: orb.append(cur);cur=apply_pi_set(pi,cur)
            d=len(orb)
            if N%d!=0: L2ok=False
            # pi^d fixes K
            pid=K
            for _ in range(d): pid=apply_pi_set(pi,pid)
            if pid!=K: L2ok=False
report(f"L1 orbit-evaporation: pi-orbit-mates share one evap time  ({count} graphs)",L1ok)
report("L2 component bundle: orbit size d | N, component is <pi^d>-invariant",L2ok)

# L5: bundle factor d^(k-1) (labeled case d=c)
A007134=[0,1,1,4,35]  # index = #vertices; connected labeled chordal
def count_single_period_c_bundles(c,k):
    # cycle type c^k (k c-cycles), n=c*k; count pi-invariant chordal graphs whose components
    # form a SINGLE pi-orbit of size exactly c (one period-c bundle)
    ct=[c]*k; pi,cyc,n=perm_of(ct); tot=0
    for adj,_ in all_invariant_graphs(ct):
        if not is_chordal(adj,n): continue
        comps=components(adj,n)
        if len(comps)!=c: continue
        # single orbit under pi?
        K=comps[0]; orb={K};cur=apply_pi_set(pi,K)
        while cur!=K: orb.add(cur);cur=apply_pi_set(pi,cur)
        if len(orb)==c and set(comps)==orb: tot+=1
    return tot
ok=True;detail=[]
for (c,k) in [(2,2),(2,3),(3,2),(4,2)]:
    got=count_single_period_c_bundles(c,k); pred=A007134[k]*(c**(k-1))
    detail.append(f"c={c},k={k}: got={got} pred=A007134({k})*{c}^{k-1}={pred}")
    if got!=pred: ok=False
report("L5 bundle factor: #(period-c bundles on c^k) = A007134(k)*c^(k-1)",ok)
for dline in detail: print("        "+dline,flush=True)

print(f"\n{'ALL LEMMA CHECKS PASS' if fails==0 else str(fails)+' FAILURES'}",flush=True)
