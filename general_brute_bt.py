#!/usr/bin/env python3
# brute oracle: fix/cfix by backtracking over orbit edges with chordality pruning
import os
_HERE=os.path.dirname(os.path.abspath(__file__))
import sys
from math import gcd
from itertools import combinations

def cycles_of(ct):
    cyc=[]; base=0
    for c in ct:
        cyc.append(list(range(base,base+c))); base+=c
    return cyc, base

def within_orbits(cyc):
    c=len(cyc); orbs=[]
    for d in range(1,c//2+1):
        # distance-d circulant edges; orbit = all of them (one orbit per distinct distance,
        # except for even c & d=c/2 the pairs are {i,i+d} with i<i+d -> c/2 distinct)
        seen=set(); orb=[]
        for i in range(c):
            j=(i+d)%c; e=(min(cyc[i],cyc[j]),max(cyc[i],cyc[j]))
            if e not in seen: seen.add(e); orb.append(e)
        orbs.append(orb)
    return orbs

def between_orbits(cycA,cycB):
    cA,cB=len(cycA),len(cycB); g=gcd(cA,cB)
    seen=[[False]*cB for _ in range(cA)]; orbs=[]
    for i in range(cA):
        for j in range(cB):
            if seen[i][j]: continue
            orb=[]; ii,jj=i,j
            while not seen[ii][jj]:
                seen[ii][jj]=True
                u,v=cycA[ii],cycB[jj]; orb.append((min(u,v),max(u,v)))
                ii=(ii+1)%cA; jj=(jj+1)%cB
            orbs.append(orb)
    return orbs

def is_chordal_set(adj, verts):
    alive=list(verts)
    aset=set(alive)
    while aset:
        found=None
        for v in aset:
            Nv=adj[v]&aset; Nv=Nv-{v}
            ok=True
            for u in Nv:
                if (adj[u]&Nv)!=(Nv-{u}): ok=False; break
            if ok: found=v; break
        if found is None: return False
        aset.discard(found)
    return True

def solve(ct):
    cyc,n=cycles_of(ct)
    adj=[set() for _ in range(n)]
    placed=[]
    fixc=[0]; cfixc=[0]
    r=len(cyc)
    def rec(k):
        if k==r:
            fixc[0]+=1
            # connectivity
            seen={0}; st=[0]
            while st:
                v=st.pop()
                for u in adj[v]:
                    if u not in seen: seen.add(u); st.append(u)
            if len(seen)==n: cfixc[0]+=1
            return
        C=cyc[k]
        neo=within_orbits(C)
        for j in range(k):
            neo+=between_orbits(C,cyc[j])
        m=len(neo)
        newverts=set(placed)|set(C)
        for mask in range(1<<m):
            added=[]
            mm=mask; bit=0
            while mm:
                if mm&1:
                    for (u,v) in neo[bit]:
                        adj[u].add(v); adj[v].add(u); added.append((u,v))
                mm>>=1; bit+=1
            if is_chordal_set(adj,newverts):
                placed.extend(C); rec(k+1)
                del placed[len(placed)-len(C):]
            for (u,v) in added: adj[u].discard(v); adj[v].discard(u)
    rec(0)
    return fixc[0],cfixc[0],n

if __name__=="__main__":
    if len(sys.argv)>1:
        ct=[int(x) for x in sys.argv[1].split(",")]
        import time; t0=time.time()
        f,c,n=solve(ct)
        print(f"lambda={ct} n={n}: fix={f} cfix={c}  [{time.time()-t0:.2f}s]")
    else:
        print("validate involution 2^m (expect cfix 1,7,154,7289):")
        for m in range(1,5):
            f,c,n=solve([2]*m); print(f"  2^{m}: fix={f} cfix={c}")
        print("validate identity 1^n (expect cfix 1,1,4,35,541):")
        for nn in range(1,6):
            f,c,n=solve([1]*nn); print(f"  1^{nn}: fix={f} cfix={c}")
        print("mixed/large-cycle spot checks:")
        for ct in [[3,2,1],[4,4],[3,3,3,3,3,1],[8,8],[16]]:
            import time; t0=time.time()
            f,c,n=solve(ct); print(f"  {ct}: fix={f} cfix={c}  [{time.time()-t0:.2f}s]")
