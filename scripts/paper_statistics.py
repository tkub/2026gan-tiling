#!/usr/bin/env python3
"""
paper_statistics.py -- Prints the statistics quoted in the paper,
recomputed from the data files in data/.

Covers: catalog-wide n_til statistics, the per-class verification of
the n_til-max rule and within-class Spearman correlations (Table III),
descriptor correlations (Sec. III E), the log(n_til) relation
(Sec. IV B), the MLIP-DFT comparison statistics for the
fixed-placement and cross-class DFT validation sets
(Sec. III C, III D), and the local frustration count n_adj
(Sec. III F): within-class correlations, energy per adjacent pair,
mediation of the n_til correlation, and the n_adj = 0 census.

Run from the repository root:  python3 scripts/paper_statistics.py
"""

import csv
import math
from collections import Counter, defaultdict


def median(v):
    s = sorted(v)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sx = math.sqrt(sum((a - mx) ** 2 for a in x))
    sy = math.sqrt(sum((b - my) ** 2 for b in y))
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / (sx * sy)


def ranks(v):
    idx = sorted(range(len(v)), key=lambda i: v[i])
    rk = [0.0] * len(v)
    i = 0
    while i < len(idx):
        j = i
        while j + 1 < len(idx) and v[idx[j + 1]] == v[idx[i]]:
            j += 1
        for k in range(i, j + 1):
            rk[idx[k]] = (i + j) / 2
        i = j + 1
    return rk


def spearman(x, y):
    return pearson(ranks(x), ranks(y))


def main():
    # -- catalog --
    hist = Counter()
    ntmax = {}
    at_max = defaultdict(int)
    with open('data/position_catalog_kG3.csv') as f:
        for row in csv.DictReader(f):
            nt = int(row['n_tilings'])
            c = int(row['ga_orbit']) + 1
            hist[nt] += 1
            if nt > ntmax.get(c, -1):
                ntmax[c], at_max[c] = nt, 1
            elif nt == ntmax[c]:
                at_max[c] += 1
    total = sum(hist.values())
    print("== catalog ==")
    print(f"configurations: {total:,}")
    print(f"n_til = 1: {hist[1]:,} ({100 * hist[1] / total:.1f}%), "
          f"mean {sum(k * v for k, v in hist.items()) / total:.2f}, "
          f"max {max(hist)}")
    print(f"n_til-max candidates: {sum(at_max.values())} "
          f"(per class: {dict(sorted(at_max.items()))})")

    # -- relaxed sample --
    rel = []
    with open('data/descriptors_relaxed.csv') as f:
        for r in csv.DictReader(f):
            rel.append((int(r['config_id']), int(r['ga_orbit']) + 1,
                        int(r['n_til']), float(r['energy_eV']),
                        float(r['I_sh']), float(r['I_L2']),
                        float(r['V_HH']), float(r['d_GaH_mean'])))
    E0 = min(r[3] for r in rel)
    dE = [(r[3] - E0) * 1000 for r in rel]
    nt = [r[2] for r in rel]
    print(f"\n== relaxed sample (n = {len(rel)}) ==")
    print(f"global Spearman rho(n_til, E) = {spearman(nt, dE):.2f}")
    print(f"Pearson r(n_til, E) = {pearson(nt, dE):.2f}, "
          f"r(ln n_til, E) = {pearson([math.log(v) for v in nt], dE):.2f}")

    byc = defaultdict(list)
    for r, de in zip(rel, dE):
        byc[r[1]].append((r[2], de, r[0]))
    print("\nclass  n    ntmax  E_min config  n_til(E_min)  gap(meV)  rho")
    rhos, slopes = [], []
    for c in sorted(byc):
        lst = byc[c]
        x = [a[0] for a in lst]
        y = [a[1] for a in lst]
        rho = spearman(x, y)
        rhos.append(rho)
        emin = min(lst, key=lambda a: a[1])
        bam = min((a for a in lst if a[0] == ntmax[c]),
                  key=lambda a: a[1])
        lx = [math.log(v) for v in x]
        mx, my = sum(lx) / len(lx), sum(y) / len(y)
        sxx = sum((v - mx) ** 2 for v in lx)
        if sxx > 1e-9:
            slopes.append(sum((a - mx) * (b - my)
                              for a, b in zip(lx, y)) / sxx)
        print(f"{c:3d} {len(lst):5d} {ntmax[c]:5d} {emin[2]:12d} "
              f"{emin[0]:8d} {bam[1] - emin[1]:11.1f} {rho:+9.2f}")
    rhos.sort()
    slopes.sort()
    print(f"within-class rho: median {median(rhos):.2f}, "
          f"range {rhos[0]:.2f} .. {rhos[-1]:.2f}")
    # robustness against nonuniform sampling density: collapse each
    # class to one point per n_til value (median energy per value)
    rb = []
    for c in sorted(byc):
        bins = defaultdict(list)
        for ntv, de, _ in byc[c]:
            bins[ntv].append(de)
        bx = sorted(bins)
        if len(bx) > 2:
            bmed = [median(bins[v]) for v in bx]
            rb.append(spearman(bx, bmed))
    rb.sort()
    print(f"within-class rho, one point per n_til value: "
          f"median {median(rb):.2f}, "
          f"range {rb[0]:.2f} .. {rb[-1]:.2f}")
    print(f"slope of E on ln(n_til): median {median(slopes):.0f} "
          f"meV per e-fold, range {slopes[0]:.0f} .. {slopes[-1]:.0f}")
    # robustness margin: E(best config with n_til < class max) minus
    # E(best n_til-max config), per class
    print("\nmargin of the n_til-max configuration (meV):")
    for c in sorted(byc):
        lst = byc[c]
        best_max = min(de for ntv, de, _ in lst if ntv == ntmax[c])
        best_rest = min(de for ntv, de, _ in lst if ntv < ntmax[c])
        print(f"  class {c:2d}: {best_rest - best_max:+7.1f}")

    # -- descriptor correlations (Sec. III E) --
    print("\n== geometric descriptors vs energy ==")
    for name, idx in (('I_L2 (I_iso)', 5), ('I_sh', 4), ('V_HH', 6),
                      ('d_GaH_mean', 7)):
        print(f"r({name}, E) = {pearson([r[idx] for r in rel], dE):+.2f}")

    # -- fixed-placement DFT validation set --
    rows = []
    with open('data/dft_uma_mapped.csv') as f:
        for r in csv.DictReader(f):
            rows.append((float(r['dE_DFT_eV']), float(r['dE_UMA_eV']),
                         int(r['n_til'])))
    dd = [r[0] for r in rows]
    du = [r[1] for r in rows]
    print(f"\n== fixed-placement DFT validation set (n = {len(rows)}) ==")
    print(f"Spearman rho = {spearman(du, dd):.3f}, "
          f"Pearson r = {pearson(du, dd):.3f}, "
          f"MAE = {sum(abs(a - b) for a, b in zip(du, dd)) / len(dd) * 1000:.0f} meV")
    med = median
    for k in (0, 1, 2, 4):
        sub = [r[0] for r in rows if r[2] == k]
        if sub:
            print(f"n_til = {k}: n = {len(sub):3d}, "
                  f"median dE_DFT = {med(sub):.3f} eV")
    # linear calibration: DFT = a * UMA + b
    n = len(du)
    mx, my = sum(du) / n, sum(dd) / n
    a = (sum((u - mx) * (d - my) for u, d in zip(du, dd))
         / sum((u - mx) ** 2 for u in du))
    b = my - a * mx
    res = [d - (a * u + b) for u, d in zip(du, dd)]
    print(f"linear calibration: dE_DFT = {a:.3f} * dE_UMA "
          f"+ {b * 1000:.0f} meV; "
          f"calibrated MAE = {sum(abs(e) for e in res) / n * 1000:.0f} meV")

    # -- cross-class DFT validation set --
    dft = {}
    with open('data/dft_results.csv') as f:
        for r in csv.DictReader(f):
            dft[int(r['config_id'])] = float(r['E_dft_total_eV'])
    mlip = {r[0]: r[3] for r in rel}
    common = sorted(set(dft) & set(mlip))
    x = [mlip[c] for c in common]
    y = [dft[c] for c in common]
    n = len(common)
    mx, my = sum(x) / n, sum(y) / n
    resid = [(a - mx) - (b - my) for a, b in zip(x, y)]
    print(f"\n== cross-class DFT validation set (n = {n}) ==")
    print(f"Spearman rho = {spearman(x, y):.3f}, "
          f"Pearson r = {pearson(x, y):.3f}, "
          f"MAE = {sum(abs(d) for d in resid) / n * 1000:.0f} meV")
    dmin = min(y)
    low = [(mlip[c], dft[c]) for c in common
           if (dft[c] - dmin) * 1000 < 200]
    lmx = sum(a for a, _ in low) / len(low)
    lmy = sum(b for _, b in low) / len(low)
    lmae = sum(abs((a - lmx) - (b - lmy))
               for a, b in low) / len(low) * 1000
    print(f"low-energy subset (dE_DFT < 200 meV): n = {len(low)}, "
          f"MAE = {lmae:.0f} meV")

    # -- n_til-max candidate verification by DFT (Secs. III D, IV D) --
    cls_d, nt_d = {}, {}
    with open('data/position_catalog_kG3.csv') as f:
        for row in csv.DictReader(f):
            cid = int(row['config_id'])
            if cid in dft:
                cls_d[cid] = int(row['ga_orbit']) + 1
                nt_d[cid] = int(row['n_tilings'])
    cand = {cid for cid in dft if nt_d[cid] == ntmax[cls_d[cid]]}
    assert len(cand) == 24 and set(dft) - cand == {21900}
    byc_d = defaultdict(list)
    for cid in dft:
        byc_d[cls_d[cid]].append(cid)
    mins = {c: min(l, key=lambda q: dft[q]) for c, l in byc_d.items()}
    n_ok = sum(mins[c] in cand for c in mins)
    print("\n== n_til-max candidate verification by DFT ==")
    print(f"all 24 candidates have DFT energies; the lowest DFT config "
          f"of the class is a candidate in {n_ok}/14 classes")
    assert n_ok == 13 and mins[2] == 21900
    ee = sorted(dft[q] for q in mins.values())
    print(f"spread of the 14 class minima (DFT): "
          f"{(ee[-1] - ee[0]) * 1000:.0f} meV")
    print(f"E_DFT(109117) - E_DFT(21900) = "
          f"{(dft[109117] - dft[21900]) * 1000:+.1f} meV")
    c5 = [dft[c] for c in byc_d[5]]
    print(f"Class-5 candidate span (DFT): "
          f"{(max(c5) - min(c5)) * 1000:.0f} meV")


    # -- placement-chain geometry (Sec. III F, Class-2 exception) --
    import itertools
    canon2 = {}
    with open('data/position_catalog_kG3.csv') as f:
        for row in csv.DictReader(f):
            k = int(row['ga_orbit'])
            canon2.setdefault(
                k, tuple(int(x) for x in row['ga_faces'].split(';')))
    chains = {}
    for k, faces in sorted(canon2.items()):
        pos = [((q // 12), (q % 12) // 2) for q in faces]
        sp = None
        for perm in itertools.permutations(pos):
            v = [(((perm[(i + 1) % 3][0] - perm[i][0]) % 6),
                  ((perm[(i + 1) % 3][1] - perm[i][1]) % 6))
                 for i in range(3)]
            if v[0] == v[1] == v[2]:
                dq, dr = v[0]
                sp = min((dq + a2) ** 2 + (dr + b2) ** 2
                         + (dq + a2) * (dr + b2)
                         for a2 in (-6, 0, 6) for b2 in (-6, 0, 6))
        chains[k + 1] = sp
    assert chains[2] == 4 and chains[1] == 12
    assert all(chains[c] is None for c in range(3, 15))
    print("\n== placement geometry ==")
    print("Class 2 is the only periodic Ga chain with spacing 2a; "
          "Class 1 is collinear at 2*sqrt(3)*a; Classes 3-14 are "
          "non-collinear")

    # -- local frustration count n_adj (Sec. III F) --
    print("\n== local frustration count n_adj ==")
    tri = [f for f in range(72) if f % 2 == 1]

    def _cent(f):
        t2 = f % 2
        k = f // 2
        q, r2 = k // 6, k % 6
        if t2 == 0:
            vv = [(q, r2), (q + 1, r2), (q, r2 + 1)]
        else:
            vv = [(q + 1, r2), (q + 1, r2 + 1), (q, r2 + 1)]
        s3 = math.sqrt(3) / 2
        return (sum(a2 + b2 / 2 for a2, b2 in vv) / 3,
                -sum(b2 for a2, b2 in vv) / 3 * s3)

    def _dist(a2, b2):
        s3 = math.sqrt(3) / 2
        v1, v2 = (6.0, 0.0), (3.0, -6 * s3)
        c1, c2 = _cent(a2), _cent(b2)
        return min(math.hypot(c1[0] - c2[0] + i2 * v1[0] + j2 * v2[0],
                              c1[1] - c2[1] + i2 * v1[1] + j2 * v2[1])
                   for i2 in (-1, 0, 1) for j2 in (-1, 0, 1))

    nn = {t2: frozenset(u for u in tri if u != t2
                        and abs(_dist(t2, u) - 1.0) < 1e-3) for t2 in tri}
    gadj = {g: frozenset(t2 for t2 in tri
                         if abs(_dist(g, t2) - 1 / math.sqrt(3)) < 1e-3)
            for g in range(0, 72, 2)}

    def nadj_of(ga2, h2):
        ex = set()
        for g in ga2:
            ex |= gadj[g]
        bare = set(t2 for t2 in tri if t2 not in h2 and t2 not in ex)
        return sum(len(nn[v] & bare) for v in bare) // 2

    faces = {}
    nadj_all = {}
    ntil_all = {}
    cls_all = {}
    nadj0 = 0
    with open('data/position_catalog_kG3.csv') as f:
        for row in csv.DictReader(f):
            cid = int(row['config_id'])
            ga2 = [int(x2) for x2 in row['ga_faces'].split(';')]
            h2 = set(int(x2) for x2 in row['h_faces'].split(';'))
            na = nadj_of(ga2, h2)
            nadj_all[cid] = na
            ntil_all[cid] = int(row['n_tilings'])
            cls_all[cid] = int(row['ga_orbit']) + 1
            if na == 0:
                nadj0 += 1
    print(f"catalog: n_adj = 0 for {nadj0} of {len(nadj_all):,} configs")

    byc2 = defaultdict(list)
    for r in rel:
        byc2[r[1]].append((nadj_all[r[0]], (r[3] - E0) * 1000, r[2], r[0]))
    rhos2, slopes2 = [], []
    for c in sorted(byc2):
        v = byc2[c]
        x2 = [a2 for a2, b2, n2, c2 in v]
        y2 = [b2 for a2, b2, n2, c2 in v]
        rhos2.append(spearman(x2, y2))
        mx2, my2 = sum(x2) / len(x2), sum(y2) / len(y2)
        slopes2.append(sum((a2 - mx2) * (b2 - my2)
                           for a2, b2 in zip(x2, y2))
                       / sum((a2 - mx2) ** 2 for a2 in x2))
        emin_c = min(v, key=lambda t2: t2[1])
        assert emin_c[0] == 0, f"class {c}: minimum has n_adj != 0"
    rhos2.sort()
    slopes2.sort()
    print(f"within-class rho(n_adj, E): median "
          f"{median(rhos2):+.2f}, "
          f"range {rhos2[0]:+.2f} .. {rhos2[-1]:+.2f}")
    print(f"within-class slope: median {median(slopes2):.0f} "
          f"meV per adjacent pair, "
          f"range {slopes2[0]:.0f} .. {slopes2[-1]:.0f}")
    print("every class minimum has n_adj = 0: OK")
    # mediation
    cells1 = defaultdict(list)
    cells2 = defaultdict(list)
    for c in byc2:
        for a2, b2, n2, cid in byc2[c]:
            cells1[(c, n2)].append((a2, b2))
            cells2[(c, a2)].append((n2, b2))

    def cellmed(cells):
        rr = []
        for k2, v2 in cells.items():
            if len(v2) < 6:
                continue
            xv = [a2 for a2, b2 in v2]
            if len(set(xv)) < 2:
                continue  # skip cells where the feature is constant
            r2 = spearman(xv, [b2 for a2, b2 in v2])
            if r2 == r2:
                rr.append(r2)
        rr.sort()
        return len(rr), median(rr)

    n1c, m1c = cellmed(cells1)
    n2c, m2c = cellmed(cells2)
    print(f"rho(n_adj, E | class, n_til): median {m1c:+.2f} "
          f"({n1c} cells)")
    print(f"rho(n_til, E | class, n_adj): median {m2c:+.2f} "
          f"({n2c} cells)")
    # DFT (fixed-placement, all 685 incl. EC-incompatible)
    xs, ys = [], []
    with open('data/dft_uma_mapped.csv') as f:
        for row2 in csv.DictReader(f):
            ga2 = [int(x2) for x2 in row2['ga_faces'].split(';')]
            h2 = set(int(x2) for x2 in row2['h_faces'].split(';'))
            xs.append(nadj_of(ga2, h2))
            ys.append(float(row2['dE_DFT_eV']))
    mx2, my2 = sum(xs) / len(xs), sum(ys) / len(ys)
    sl2 = (sum((a2 - mx2) * (b2 - my2) for a2, b2 in zip(xs, ys))
           / sum((a2 - mx2) ** 2 for a2 in xs))
    print(f"DFT (n = {len(xs)}): rho(n_adj, dE_DFT) = "
          f"{spearman(xs, ys):+.2f}, slope = {1000 * sl2:.0f} meV/pair")

    # -- complete n_adj = 0 set and untileable arrangements (Sec. III F) --
    print("\n== complete n_adj = 0 set (102) and untileable arrangements ==")
    res167 = {int(r['config_id']): float(r['energy_eV'])
              for r in csv.DictReader(open('data/relax_all167_out.csv'))
              if r['converged'] == '1'}
    census = []
    clsmin = {}
    for r0 in rel:
        c0 = r0[1]
        clsmin[c0] = min(clsmin.get(c0, 1e9), r0[3])
        if nadj_all.get(r0[0]) == 0:
            census.append((c0, ntil_all[r0[0]], r0[3]))
    for r in csv.DictReader(open('data/relax_needed_nadj0.csv')):
        cid = int(r['config_id'])
        c0 = int(r['ga_orbit']) + 1
        e0 = res167[cid]
        census.append((c0, int(r['n_tilings']), e0))
        assert e0 >= clsmin[c0] - 1e-9
    assert len(census) == 102
    byc0 = defaultdict(list)
    for c0, nt0, e0 in census:
        byc0[c0].append((nt0, e0))
    rr, hits, gaps = [], 0, []
    for c0, v0 in byc0.items():
        rr.append(spearman([a0 for a0, b0 in v0],
                           [b0 for a0, b0 in v0]))
        ntm = max(a0 for a0, b0 in v0)
        bam = min(b0 for a0, b0 in v0 if a0 == ntm)
        if abs(bam - min(b0 for a0, b0 in v0)) < 1e-9:
            hits += 1
            rest = min((b0 for a0, b0 in v0 if a0 < ntm), default=None)
            if rest is not None:
                gaps.append(1000 * (rest - bam))
    rr.sort()
    gaps.sort()
    print(f"within-class rho(n_til, E) on the n_adj = 0 set: median "
          f"{median(rr):+.2f} (all negative: {all(x < 0 for x in rr)})")
    print(f"n_til-max = minimum of the set: {hits}/14; margins "
          f"{gaps[0]:.0f}-{gaps[-1]:.0f} meV")
    esc = {900037, 900040, 900041, 900068, 900084}
    marg = []
    for r in csv.DictReader(open('data/relax_nadj0_untileable.csv')):
        cid = int(r['config_id'])
        if cid in esc:
            continue
        c0 = int(r['ga_orbit']) + 1
        m0 = 1000 * (res167[cid] - min(
            clsmin[c0],
            min((res167[int(r2['config_id'])]
                 for r2 in csv.DictReader(open('data/relax_needed_nadj0.csv'))
                 if int(r2['ga_orbit']) + 1 == c0), default=1e9)))
        marg.append(m0)
    marg.sort()
    print(f"untileable (83 well-behaved): +{marg[0]:.0f}..+{marg[-1]:.0f} "
          f"meV above class minima (median +{median(marg):.0f}); "
          f"5 escaped structures excluded")
    # combined margins quoted in Sec. III D
    allrel = defaultdict(list)
    for r0 in rel:
        allrel[r0[1]].append((r0[2], r0[3]))
    for r in csv.DictReader(open('data/relax_needed_nadj0.csv')):
        allrel[int(r['ga_orbit']) + 1].append(
            (int(r['n_tilings']), res167[int(r['config_id'])]))
    cm = []
    for c0, v0 in allrel.items():
        ntm = max(a0 for a0, b0 in v0)
        m0 = 1000 * (min(b0 for a0, b0 in v0 if a0 < ntm)
                     - min(b0 for a0, b0 in v0 if a0 == ntm))
        if m0 > 0:
            cm.append(m0)
    cm.sort()
    print(f"rule margins over all relaxed configs: {cm[0]:.0f}-{cm[-1]:.0f} "
          f"meV, median {median(cm):.0f} (13 classes)")


if __name__ == '__main__':
    main()
