#!/usr/bin/env python3
"""Verify the bare-site arrangement counts of Sec. III F.

Checks, by direct enumeration of independent 9-subsets at the
canonical representative placement of each of the 14 classes:
  * the 14 placements admit 246 bare-site arrangements with n_adj = 0;
  * 88 of them (36%) admit no compatible tiling (absent from the
    catalog);
  * the remaining 158 tileable arrangements reduce, under the
    placement stabilizers, to the 102 symmetry-distinct
    configurations of the n_adj = 0 census.

Group operations are taken from enumeration.py.

Run from the repository root:
    python3 scripts/verify_nadj0_arrangements.py
"""
import csv
import math
import sys

sys.path.insert(0, 'scripts')
import enumeration as en  # noqa: E402

S3 = math.sqrt(3) / 2
V1, V2 = (6.0, 0.0), (3.0, -6 * S3)
TRI = [f for f in range(72) if f % 2 == 1]


def cent(f):
    q, r, t = f // 12, (f % 12) // 2, f % 2
    vv = ([(q, r), (q + 1, r), (q, r + 1)] if t == 0
          else [(q + 1, r), (q + 1, r + 1), (q, r + 1)])
    return (sum(a + b / 2 for a, b in vv) / 3,
            -sum(b for a, b in vv) / 3 * S3)


def dist(a, b):
    c1, c2 = cent(a), cent(b)
    return min(math.hypot(c1[0] - c2[0] + i * V1[0] + j * V2[0],
                          c1[1] - c2[1] + i * V1[1] + j * V2[1])
               for i in (-1, 0, 1) for j in (-1, 0, 1))


def main():
    _, Gp = en.build_symmetry_perms()
    NN = {t: {u for u in TRI if u != t and abs(dist(t, u) - 1.0) < 1e-3}
          for t in TRI}
    canon = {}
    rows = []
    with open('data/position_catalog_kG3.csv') as fp:
        for r in csv.DictReader(fp):
            k = int(r['ga_orbit'])
            ga = frozenset(int(x) for x in r['ga_faces'].split(';'))
            canon.setdefault(k, ga)
            assert ga == canon[k], "catalog not at canonical placements?"
            rows.append((k, frozenset(int(x)
                                      for x in r['h_faces'].split(';'))))

    tot_raw = tot_until = tot_til = tot_census = 0
    for k in sorted(canon):
        ga = canon[k]
        excl = {t for g in ga for t in TRI
                if abs(dist(g, t) - 1 / math.sqrt(3)) < 1e-3}
        allowed = [t for t in TRI if t not in excl]
        stab = [p for p in Gp
                if frozenset(p[f] for f in ga) == ga]
        # census members of this class (n_adj = 0 catalog entries)
        ch = set()
        for kk, h in rows:
            if kk != k:
                continue
            bare = [t for t in allowed if t not in h]
            nn = sum(1 for i, a in enumerate(bare)
                     for b in bare[i + 1:] if b in NN[a])
            if nn == 0:
                ch.add(h)
        # all independent 9-subsets of the 27 allowed sites
        res = []

        def dfs(idx, chosen):
            if len(chosen) == 9:
                res.append(frozenset(chosen))
                return
            if idx >= len(allowed) or len(allowed) - idx < 9 - len(chosen):
                return
            t = allowed[idx]
            if not (NN[t] & set(chosen)):
                chosen.append(t)
                dfs(idx + 1, chosen)
                chosen.pop()
            dfs(idx + 1, chosen)

        dfs(0, [])
        til = until = 0
        for bare in res:
            h = frozenset(t for t in allowed if t not in bare)
            if any(frozenset(p[f] for f in h) in ch for p in stab):
                til += 1
            else:
                until += 1
        print(f"class {k + 1:2d}: raw = {len(res):3d}, "
              f"tileable = {til:3d}, untileable = {until:3d}, "
              f"census = {len(ch):3d}, |stab| = {len(stab)}")
        tot_raw += len(res)
        tot_until += until
        tot_til += til
        tot_census += len(ch)

    print(f"TOTAL: raw = {tot_raw}, tileable = {tot_til}, "
          f"untileable = {tot_until}, census = {tot_census}")
    assert tot_raw == 246 and tot_until == 88
    assert tot_til == 158 and tot_census == 102
    print("OK: 246 arrangements = 88 untileable + 158 tileable; "
          "the 158 reduce to the 102 census configurations.")


if __name__ == '__main__':
    main()
