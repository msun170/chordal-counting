#!/usr/bin/env python3
r"""Independent cross-check: chordal graphs are closed under disjoint union, so
A048193(all) = Euler transform of A048192(connected). Verify this holds for ALL n=1..16,
INCLUDING the new a(16) -- this validates both new terms using only the connected sequence."""
from fractions import Fraction

# my Burnside-computed values n=1..20
CONN=[1,1,2,5,15,58,272,1614,11911,109539,1247691,17566431,305310547,6558690953,174688164414,5796153514484,241003010628949,12642592677074970,842762851699294393,71916937400532750123]
ALL =[1,2,4,10,27,94,393,2119,14524,126758,1392387,19109099,326005775,6905776799,181945055235,5985406996403,247178491630853,12895963060540295,855912598965399807,72786012927793961715]

def euler_transform(c):
    # c[0]=c_1,... connected counts. Returns a[1..N] all counts via
    # 1+sum a_n x^n = prod_k (1-x^k)^(-c_k). Standard integer recurrence.
    N=len(c)
    # b_n = sum_{d|n} d*c_d
    b=[0]*(N+1)
    for k in range(1,N+1):
        for m in range(k,N+1,k):
            b[m]+=k*c[k-1]
    a=[0]*(N+1); a[0]=1
    for n in range(1,N+1):
        s=b[n]+sum(b[k]*a[n-k] for k in range(1,n))
        assert s%n==0, f"non-integer at n={n}"
        a[n]=s//n
    return a[1:]

et=euler_transform(CONN)
print(f"{'n':>3} {'Euler(conn)':>16} {'A048193(all)':>16}  check")
ok=True
for i in range(len(ALL)):
    n=i+1; e=et[i]; a=ALL[i]
    good = (e==a)
    ok &= good
    print(f"{n:>3} {e:>16} {a:>16}  {'OK' if good else '*** MISMATCH ***'}")
print("\nEULER TRANSFORM CONFIRMS ALL TERMS incl a(20)!!" if ok else "\nMISMATCH - investigate")
