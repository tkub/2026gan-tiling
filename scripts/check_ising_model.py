#!/usr/bin/env python3
"""Cross-check against the data-driven Ising model of Kawka et al.

Reconstructs the Ising model of K. Kawka et al., J. Appl. Phys. 135,
225302 (2024) and evaluates it over all EC-compatible configurations
at the Class-7 Ga placement (14,896 catalog entries).  Expected
result: the unique Ising minimum is configuration 109117, the global
MLIP minimum of the paper.

Parameter values follow the reference implementation (PyAPX,
https://github.com/a-ksb/PyAPX, examples/H_GaN0001_6x6/
my_evaluator.py): J = -0.021863 eV with each interaction pair counted
twice, i.e. J_pair = -0.043726 eV with each pair counted once as
done here; h = -0.097453 eV; c = 0.818264 eV.  Interaction pairs are
restricted to the 27 allowed sites (the 9 sites adjacent to a Ga
adatom are excluded), as in the reference implementation.

Run from the repository root:
    python3 scripts/check_ising_model.py
"""
import csv
import math

J, H, C = -0.043726, -0.097453, 0.818264
S3 = math.sqrt(3) / 2


def from_id(f):
    t = f % 2
    k = f // 2
    return k // 6, k % 6, t


def xy(q, r):
    return (q + 0.5 * r, -r * S3)


def cent(f):
    q, r, t = from_id(f)
    if t == 0:
        vv = [(q, r), (q + 1, r), (q, r + 1)]
    else:
        vv = [(q + 1, r), (q + 1, r + 1), (q, r + 1)]
    return (sum(xy(a, b)[0] for a, b in vv) / 3,
            sum(xy(a, b)[1] for a, b in vv) / 3)


V1, V2 = xy(6, 0), xy(0, 6)


def dist(f, g):
    c1, c2 = cent(f), cent(g)
    return min(math.hypot(c1[0] - c2[0] + i * V1[0] + j * V2[0],
                          c1[1] - c2[1] + i * V1[1] + j * V2[1])
               for i in (-1, 0, 1) for j in (-1, 0, 1))


def main():
    tri = [f for f in range(72) if f % 2 == 1]
    row = None
    with open('data/position_catalog_kG3.csv') as fp:
        for r in csv.DictReader(fp):
            if int(r['config_id']) == 109117:
                row = r
                break
    ga = [int(x) for x in row['ga_faces'].split(';')]
    excl = {s for g in ga for s in tri if abs(dist(g, s) - 0.5774) < 1e-3}
    second = [(g, s) for g in ga for s in tri
              if abs(dist(g, s) - 1.1547) < 1e-3]
    sites = [s for s in tri if s not in excl]
    assert len(excl) == 9 and len(sites) == 27 and len(second) == 9
    nn_sites = [(a, b) for i, a in enumerate(sites) for b in sites[i + 1:]
                if abs(dist(a, b) - 1.0) < 1e-3]

    def ising(hset):
        s = {x: (1 if x in hset else -1) for x in tri}
        return (-J * sum(s[a] * s[b] for a, b in nn_sites)
                - H * sum(s[x] for _, x in second) + C)

    cand = []
    with open('data/position_catalog_kG3.csv') as fp:
        for r in csv.DictReader(fp):
            if r['ga_faces'] == row['ga_faces']:
                cand.append((int(r['config_id']),
                             frozenset(int(x)
                                       for x in r['h_faces'].split(';')),
                             int(r['n_tilings'])))
    print(f"EC-compatible configurations at the Class-7 placement: "
          f"{len(cand)}")
    res = sorted((ising(hs), cid, nt) for cid, hs, nt in cand)
    print("top 3 by Ising energy:")
    for e, cid, nt in res[:3]:
        print(f"  E_Ising = {e:+.4f} eV  config {cid}  n_til = {nt}")
    assert res[0][1] == 109117, "Ising minimum is not 109117"
    assert res[1][0] - res[0][0] > 1e-6, "minimum is not unique"
    print(f"OK: unique Ising minimum = 109117 "
          f"(gap to 2nd: {1000 * (res[1][0] - res[0][0]):.1f} meV)")


if __name__ == '__main__':
    main()
