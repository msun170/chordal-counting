# Counting Unlabeled Chordal Graphs

This repo counts chordal graphs on n nodes, up to isomorphism. It extends two OEIS sequences that
were stuck at n=15 since 2021:

- **A048193**: all chordal graphs on n nodes
- **A048192**: connected chordal graphs on n nodes

The new values for n=16 through 20 are in `b048193.txt` and `b048192.txt`. n=16-18 have been
submitted to the OEIS (A048193 is in review, A048192 follows once a draft slot opens up), and
n=19 and 20 were computed afterward.

## the idea

Counting graphs up to isomorphism comes down to Burnside's lemma: average, over all permutations of
the vertices, how many graphs each permutation keeps fixed. So the real job is counting chordal
graphs that are invariant under a given permutation. We do that with a dynamic program built on top
of the labeled chordal counting DP of Hebert-Johnson, Lokshtanov and Vigoda (arXiv:2308.09703).

The new part is a divisor-bundle trick. When a connected piece spans a cyclic orbit of size c, it
can show up as d copies for any divisor d of c, and each copy lives in a smaller symmetry world that
the same program handles recursively. The full writeup with proofs is in `paper/paper.pdf`.

## files

- `rec_dp.py` master DP, gives fix and cfix for any cycle type
- `assemble.py` the labeled HLV DP, used as the base case
- `general_brute_bt.py` brute force oracle, used to check everything
- `burnside.py` the Burnside sum, checks against the known values
- `gen_bfile.py` computes a(n) and writes the b-files
- `verify_all.py` runs every validation layer
- `lemma_checks.py` checks the structural lemmas on actual graphs
- `euler_check.py` checks that A048193 is the Euler transform of A048192
- `valtable.py` per-cycle-type checks against brute force

## run it

needs python 3, no extra libraries. every command tells you if it passed.

**`python3 lemma_checks.py`** (fast) checks the structural lemmas on actual graphs. passed if the
last line is:
```
ALL LEMMA CHECKS PASS
```

**`python3 euler_check.py`** (fast) checks A048193 = euler transform of A048192 through n=20. passed
if the last line is:
```
EULER TRANSFORM CONFIRMS ALL TERMS incl a(20)!!
```

**`python3 valtable.py`** (fast) compares fix and cfix to brute force for 23 cycle types. passed if
the last line is:
```
% ALL MATCH
```

**`python3 verify_all.py`** (about 15 minutes) the full validation, six layers. each layer prints
`[PASS]` or `[FAIL]` as it finishes; layer 5 reproduces all 15 known terms of both sequences and
layer 6 computes a(16). passed if the last line is:
```
RESULT: ALL LAYERS PASS
```

**`python3 gen_bfile.py N`** recomputes a(1..N) and writes the b-files. each n<=15 row is tagged
`known: OK`, each n>15 row is tagged `NEW`. passed if the last line is:
```
RESULT: all known terms (n<=15) reproduced
```
use `gen_bfile.py 15` for a quick full check of the known terms (a few minutes), or `gen_bfile.py 20`
to regenerate the new terms (slow, several hours; n=19 and n=20 are a few hours each on one core).
the five new rows should read exactly:
```
n=16  A048193(all)=5985406996403  A048192(conn)=5796153514484  NEW  [...]
n=17  A048193(all)=247178491630853  A048192(conn)=241003010628949  NEW  [...]
n=18  A048193(all)=12895963060540295  A048192(conn)=12642592677074970  NEW  [...]
n=19  A048193(all)=855912598965399807  A048192(conn)=842762851699294393  NEW  [...]
n=20  A048193(all)=72786012927793961715  A048192(conn)=71916937400532750123  NEW  [...]
```

## the numbers

| n | A048193 (all) | A048192 (connected) |
|---|---|---|
| 16 | 5985406996403 | 5796153514484 |
| 17 | 247178491630853 | 241003010628949 |
| 18 | 12895963060540295 | 12642592677074970 |
| 19 | 855912598965399807 | 842762851699294393 |
| 20 | 72786012927793961715 | 71916937400532750123 |

n=1 through 15 match the known OEIS values exactly.

## checks

everything is validated four independent ways:

- reproduces all 15 known terms of both sequences
- matches brute force on every cycle type small enough to enumerate
- fix(1^n) gives A058862, the labeled chordal counts
- A048193 is the Euler transform of A048192

## sha256 of the b-files

```
b048192.txt  2933d6dd7c2fb4133032158c4a5230514762f762975734e9df598f106ef0ac94
b048193.txt  95eda5ffd21902003b79624f16dd8fae021e4a00dfd3f78e7611e95824f95abb
```

## license

MIT, see `LICENSE`. if you use this, a citation is in `CITATION.cff`.
