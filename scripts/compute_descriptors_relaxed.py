#!/usr/bin/env python3
"""
compute_descriptors_relaxed.py -- Recompute n_til and the geometric
descriptors for the relaxed sample.

For each configuration in data/relax_results.csv this script computes
  n_til   (recomputed and verified against the catalog)
  I_sh    (Shannon variant of the tiling orientational isotropy)
  I_L2    (I_iso of the paper)
  V_HH    (sum of inverse H-H distances, minimum-image convention)
  d_GaH   (mean Ga-H distance)
and writes data/descriptors_relaxed.csv.

Run from the repository root:
  python3 scripts/compute_descriptors_relaxed.py
"""

import csv
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from enumeration import (  # noqa: E402
    build_pieces, enumerate_tilings, match_catalog, piece_offsets,
    ga_site_of, from_id, face_type, N, DIRS,
)

SQRT3_2 = math.sqrt(3) / 2


def face_centroid_cart(fid):
    q, r, t = from_id(fid)
    if t == 0:
        vv = [(q, r), (q + 1, r), (q, r + 1)]
    else:
        vv = [(q + 1, r), (q + 1, r + 1), (q, r + 1)]
    cx = sum(a + b / 2.0 for a, b in vv) / 3.0
    cy = sum(b * SQRT3_2 for a, b in vv) / 3.0
    return cx, cy


def distance_matrix():
    a1 = (float(N), 0.0)
    a2 = (N * 0.5, N * SQRT3_2)
    cs = [face_centroid_cart(f) for f in range(72)]
    dist = [[0.0] * 72 for _ in range(72)]
    for i in range(72):
        for j in range(i + 1, 72):
            best = float('inf')
            for n1 in (-1, 0, 1):
                for n2 in (-1, 0, 1):
                    dx = cs[j][0] - cs[i][0] + n1 * a1[0] + n2 * a2[0]
                    dy = cs[j][1] - cs[i][1] + n1 * a1[1] + n2 * a2[1]
                    best = min(best, dx * dx + dy * dy)
            dist[i][j] = dist[j][i] = math.sqrt(best)
    return dist


def main():
    print("building tilings...")
    pieces = build_pieces()
    tilings = enumerate_tilings(pieces)
    assert len(tilings) == 456

    # blocks as bitmasks + a prefilter mask of Ga-site unions
    ga_site_map = {fs: ga_site_of(fs) for fs in pieces}
    tiling_blocks = []
    tiling_gs_union = []
    for t in tilings:
        blocks = []
        gs_union = 0
        for fs in t:
            gs_bit = 1 << ga_site_map[fs]
            gs_union |= gs_bit
            tri_mask = 0
            for f in fs:
                if face_type(f) == 1:
                    tri_mask |= 1 << f
            blocks.append((gs_bit, tri_mask))
        tiling_blocks.append(blocks)
        tiling_gs_union.append(gs_union)
    # orientation counts (internal names A/B/C)
    tiling_dirs = []
    for t in tilings:
        c = {'A': 0, 'B': 0, 'C': 0}
        for fs in t:
            c[pieces[fs][0]] += 1
        tiling_dirs.append((c['A'], c['B'], c['C']))

    dist = distance_matrix()

    # for n_adj: nearest-neighbor pairs of triangle faces and the
    # triangle faces adjacent to each downward face (Ga-adatom site)
    tri_all = [f for f in range(72) if face_type(f) == 1]
    nn_tri = {t: frozenset(u for u in tri_all if u != t
                           and abs(dist[t][u] - 1.0) < 1e-3)
              for t in tri_all}
    ga_adj = {g: frozenset(t for t in tri_all
                           if abs(dist[g][t] - 1 / math.sqrt(3)) < 1e-3)
              for g in range(72) if face_type(g) == 0}

    # target configurations
    relax = {}
    with open('data/relax_results.csv') as f:
        for row in csv.DictReader(f):
            if row['energy_eV'] == 'ERROR' or row['converged'] != '1':
                continue
            relax[int(row['config_id'])] = (float(row['energy_eV']),
                                            int(row['n_tilings']))
    cat = {}
    with open('data/position_catalog_kG3.csv') as f:
        for row in csv.DictReader(f):
            cid = int(row['config_id'])
            if cid in relax:
                cat[cid] = (int(row['ga_orbit']),
                            int(row['n_tilings']),
                            tuple(int(x) for x in row['ga_faces'].split(';')),
                            tuple(int(x) for x in row['h_faces'].split(';')))
    print(f"configs: {len(relax)} relaxed, {len(cat)} matched in catalog")
    assert len(cat) == len(relax)

    out = []
    mismatch = 0
    for i, (cid, (orbit, ntil_cat, ga, h)) in enumerate(sorted(cat.items())):
        ga_mask = 0
        for g in ga:
            ga_mask |= 1 << g
        h_mask = 0
        for x in h:
            h_mask |= 1 << x
        total = [0, 0, 0]
        n_compat = 0
        for ti in range(456):
            # necessary condition: all three Ga sit on block Ga sites
            if ga_mask & ~tiling_gs_union[ti]:
                continue
            ok = True
            n_ga = 0
            for gs_bit, tri_mask in tiling_blocks[ti]:
                if gs_bit & ga_mask:
                    if tri_mask & h_mask:
                        ok = False
                        break
                    n_ga += 1
                else:
                    if bin(tri_mask & h_mask).count('1') != 3:
                        ok = False
                        break
            if ok and n_ga == 3:
                n_compat += 1
                dirs = tiling_dirs[ti]
                for k in range(3):
                    total[k] += dirs[k]
        if n_compat != ntil_cat:
            mismatch += 1
        s = sum(total)
        fbar = [x / s for x in total] if s else [1 / 3] * 3
        ish = 0.0
        for p in fbar:
            if p > 0:
                ish -= p * math.log(p)
        ish /= math.log(3)
        il2 = 1.0 - 1.5 * sum((p - 1 / 3) ** 2 for p in fbar)
        vhh = 0.0
        for a in range(len(h)):
            for b in range(a + 1, len(h)):
                vhh += 1.0 / dist[h[a]][h[b]]
        dgah = sum(dist[g][x] for g in ga for x in h) / (len(ga) * len(h))
        # n_adj: number of nearest-neighbor pairs of bare sites
        # (bare site = triangle face with no H, not adjacent to a Ga adatom)
        excl = set()
        for g in ga:
            excl |= ga_adj[g]
        bare = set(t for t in tri_all if t not in h and t not in excl)
        assert len(bare) == 9
        nadj = sum(len(nn_tri[v] & bare) for v in bare) // 2
        out.append((cid, orbit, ntil_cat, relax[cid][0],
                    ish, il2, vhh, dgah, nadj))
        if (i + 1) % 500 == 0:
            print(f"  {i + 1}/{len(cat)}")

    print(f"n_til mismatches vs catalog: {mismatch}")
    assert mismatch == 0

    with open('data/descriptors_relaxed.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['config_id', 'ga_orbit', 'n_til', 'energy_eV',
                    'I_sh', 'I_L2', 'V_HH', 'd_GaH_mean', 'n_adj'])
        for row in out:
            w.writerow([row[0], row[1], row[2], f"{row[3]:.6f}",
                        f"{row[4]:.6f}", f"{row[5]:.6f}",
                        f"{row[6]:.4f}", f"{row[7]:.4f}", row[8]])
    print(f"wrote data/descriptors_relaxed.csv ({len(out)} rows)")


if __name__ == '__main__':
    main()
