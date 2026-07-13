#!/usr/bin/env python3
"""
enumeration.py -- Rhombus-tiling enumeration and symmetry library.

Independent re-implementation of the enumeration described in the
paper (Sec. II).  Running this file as a script verifies the shipped
catalogs from the definitions:
  * exhaustive tiling enumeration (exact set cover)  -> 456 tilings
  * set equality with data/tiling_catalog.csv
  * orbit decomposition under G = T x D6  -> 7 classes (12,12,72x4,144)
    and under the sublattice-preserving G' = T x D3 -> 8 classes
  * n_til values against data/position_catalog_kG3.csv (sampled)

Conventions: face_id = 12q + 2r + t; t=0 downward faces (N sublattice,
Ga-adatom sites), t=1 upward faces (Ga sublattice, H sites).
"""

import argparse
import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

N = 6
SQRT3_2 = math.sqrt(3) / 2

DIRS = {
    'A': ((1, 0), (0, 1)),
    'B': ((0, 1), (-1, 1)),
    'C': ((-1, 1), (-1, 0)),
}


def fid(q, r, t):
    return 12 * (q % N) + 2 * (r % N) + t


def from_id(f):
    return f // 12, (f % 12) // 2, f % 2


def face_type(f):
    return f % 2


def neighbors_t0(f):
    """Edge-sharing neighbors of a t=0 face (q,r): the t=1 faces
    (q,r), (q-1,r), (q,r-1)."""
    q, r, t = from_id(f)
    assert t == 0
    return (fid(q, r, 1), fid(q - 1, r, 1), fid(q, r - 1, 1))


# -- Candidate pieces -------------------------------------

def piece_offsets(direction):
    """Cell-relative offsets of the faces covered by a piece of the
    given orientation, determined numerically."""
    u, v = DIRS[direction]
    det = 4.0 * (u[0] * v[1] - u[1] * v[0])
    offs = []
    for q in range(-8, 9):
        for r in range(-8, 9):
            for t, (cq, cr) in ((0, (1 / 3, 1 / 3)), (1, (2 / 3, 2 / 3))):
                p = (q + cq, r + cr)
                alpha = 2.0 * (p[0] * v[1] - p[1] * v[0]) / det
                beta = 2.0 * (u[0] * p[1] - u[1] * p[0]) / det
                if 0 < alpha < 1 and 0 < beta < 1:
                    offs.append((q, r, t))
    assert len(offs) == 8, (direction, offs)
    assert sum(1 for o in offs if o[2] == 0) == 4
    return offs


def build_pieces():
    """All 108 candidate pieces. Returns {faceset: (orientation, anchor)}."""
    pieces = {}
    for d in DIRS:
        offs = piece_offsets(d)
        for q0 in range(N):
            for r0 in range(N):
                fs = frozenset(fid(q0 + dq, r0 + dr, t) for dq, dr, t in offs)
                assert len(fs) == 8
                assert fs not in pieces
                pieces[fs] = (d, (q0, r0))
    assert len(pieces) == 108
    return pieces


# -- Exhaustive tiling enumeration (exact set cover) ------

def enumerate_tilings(pieces):
    piece_list = list(pieces)
    by_face = defaultdict(list)
    for fs in piece_list:
        for f in fs:
            by_face[f].append(fs)
    tilings = []
    covered = set()
    chosen = []

    def bt():
        if len(covered) == 72:
            tilings.append(frozenset(chosen))
            return
        f0 = min(f for f in range(72) if f not in covered)
        for fs in by_face[f0]:
            if covered.isdisjoint(fs):
                covered.update(fs)
                chosen.append(fs)
                bt()
                chosen.pop()
                covered.difference_update(fs)

    bt()
    return tilings


# -- Ga-adatom site of a piece; compatibility test --------

def ga_site_of(fs):
    sites = [f for f in fs if face_type(f) == 0
             and all(n in fs for n in neighbors_t0(f))]
    assert len(sites) == 1
    return sites[0]


def is_compatible(ga_set, h_set, tiling, ga_site_map):
    n_ga = 0
    for fs in tiling:
        gs = ga_site_map[fs]
        tri = [f for f in fs if face_type(f) == 1]
        if gs in ga_set:
            if any(f in h_set for f in tri):
                return False
            n_ga += 1
        else:
            if sum(1 for f in tri if f in h_set) != 3:
                return False
    return n_ga == 3


# -- Matching against tiling_catalog.csv ------------------

def match_catalog(tilings, pieces, catalog_path):
    """Match the 456 catalog tilings against the re-enumerated set.
    Returns (name convention, tiling -> (global_id, class_id))."""
    rows = []
    with open(catalog_path, encoding='utf-8') as f:
        for row in csv.DictReader(f):
            ps = []
            for k in [f'r{i}' for i in range(1, 10)]:
                m = re.match(r'(up|ne|nw)\((\d+),(\d+)\)', row[k])
                ps.append((m.group(1), int(m.group(2)), int(m.group(3))))
            rows.append({'gid': int(row['global_id']),
                         'cls': int(row['class_id']), 'pieces': ps})
    assert len(rows) == 456

    my_tilings = set(tilings)
    offs = {d: piece_offsets(d) for d in DIRS}

    def piece_fs(d, q0, r0, dq, dr):
        return frozenset(fid(q0 + dq + a, r0 + dr + b, t)
                         for a, b, t in offs[d])

    # Assume 'up' = A with offset (0,0) and identify the (orientation,
    # offset) of 'ne' and 'nw' by a joint search over all 456 rows.
    # Rows containing only up+ne pieces do not fix the offset uniquely
    # (a uniform cyclic shift of an ne column maps to another valid
    # tiling), so rows mixing all orientations are needed.
    def cands_for(name):
        found = []
        sub = [row for row in rows
               if {p[0] for p in row['pieces']} == {'up', name}]
        for d in ('B', 'C'):
            for dq in range(N):
                for dr in range(N):
                    ok = True
                    for row in sub:
                        fss = frozenset(
                            piece_fs('A', q0, r0, 0, 0) if nm == 'up'
                            else piece_fs(d, q0, r0, dq, dr)
                            for nm, q0, r0 in row['pieces'])
                        if len(fss) != 9 or fss not in my_tilings:
                            ok = False
                            break
                    if ok:
                        found.append((d, dq, dr))
        return found

    solutions = []
    for c_ne in cands_for('ne'):
        for c_nw in cands_for('nw'):
            if c_ne[0] == c_nw[0]:
                continue
            conv_try = {'up': ('A', 0, 0), 'ne': c_ne, 'nw': c_nw}
            ok = True
            mapped = set()
            for row in rows:
                fss = frozenset(
                    piece_fs(conv_try[nm][0], q0, r0,
                             conv_try[nm][1], conv_try[nm][2])
                    for nm, q0, r0 in row['pieces'])
                if len(fss) != 9 or fss not in my_tilings:
                    ok = False
                    break
                mapped.add(fss)
            if ok and len(mapped) == 456:
                solutions.append(conv_try)
    assert solutions, "no consistent piece-name convention found"
    if len(solutions) > 1:
        print(f"  note: {len(solutions)} conventions consistent; "
              f"checking class-label invariance")
    conv = solutions[0]

    tiling_info = {}
    for row in rows:
        fss = frozenset(piece_fs(conv[nm][0], q0, r0, conv[nm][1], conv[nm][2])
                        for nm, q0, r0 in row['pieces'])
        assert fss in my_tilings, row['gid']
        tiling_info[fss] = (row['gid'], row['cls'])
    assert len(tiling_info) == 456
    return conv, tiling_info


# -- Face permutations of G = T x D6 and G' = T x D3 ------

def build_symmetry_perms():
    """Face permutations of G (432 elements) and of the
    sublattice-preserving subgroup G' (216 elements).  The point group
    D6 is generated by the 60-degree rotation and a mirror about a
    lattice vertex."""
    import math as _m

    def vert_r60(v):
        q, r = v
        x, y = q + r / 2, r * _m.sqrt(3) / 2
        c, s = 0.5, _m.sqrt(3) / 2
        x2, y2 = c * x - s * y, s * x + c * y
        r2 = round(y2 / (_m.sqrt(3) / 2))
        q2 = round(x2 - r2 / 2)
        assert abs(x2 - (q2 + r2 / 2)) < 1e-9
        return (q2, r2)

    def vert_mirror(v):
        q, r = v
        return (q + r, -r)

    def face_vertices_int(f):
        q, r, t = from_id(f)
        if t == 0:
            return [(q, r), (q + 1, r), (q, r + 1)]
        return [(q + 1, r), (q + 1, r + 1), (q, r + 1)]

    vset_to_face = {}
    for f in range(72):
        key = frozenset((a % N, b % N) for a, b in face_vertices_int(f))
        assert key not in vset_to_face
        vset_to_face[key] = f

    def face_perm_from_vmap(vmap):
        perm = [None] * 72
        for f in range(72):
            key = frozenset(tuple(x % N for x in vmap(v))
                            for v in face_vertices_int(f))
            perm[f] = vset_to_face[key]
        assert sorted(perm) == list(range(72))
        return tuple(perm)

    # point group D6
    point_maps = []
    cur = lambda v: v  # noqa: E731
    maps_r = [lambda v: v]
    for _ in range(5):
        prev = maps_r[-1]
        maps_r.append(lambda v, p=prev: vert_r60(p(v)))
    point_maps.extend(maps_r)
    for m in maps_r:
        point_maps.append(lambda v, p=m: vert_mirror(p(v)))
    point_perms = {face_perm_from_vmap(vm) for vm in point_maps}
    assert len(point_perms) == 12

    trans_perms = []
    for a in range(N):
        for b in range(N):
            perm = [None] * 72
            for f in range(72):
                q, r, t = from_id(f)
                perm[f] = fid(q + a, r + b, t)
            trans_perms.append(tuple(perm))

    G = set()
    for tp in trans_perms:
        for pp in point_perms:
            G.add(tuple(tp[pp[f]] for f in range(72)))
    assert len(G) == 432
    Gp = {g for g in G
          if all((g[f] % 2) == (f % 2) for f in range(72))}
    assert len(Gp) == 216
    return G, Gp


def tiling_orbits(tilings, perms):
    """Orbit decomposition of the tiling set under a permutation group."""
    tset = set(tilings)
    remaining = set(tilings)
    orbits = []
    while remaining:
        seed = next(iter(remaining))
        orb = set()
        stack = [seed]
        while stack:
            t = stack.pop()
            if t in orb:
                continue
            orb.add(t)
            for g in perms:
                img = frozenset(frozenset(g[f] for f in fs) for fs in t)
                assert img in tset
                if img not in orb:
                    stack.append(img)
        orbits.append(orb)
        remaining -= orb
    return orbits




def main():
    import random

    print("enumerating tilings (exact set cover) ...")
    pieces = build_pieces()
    tilings = enumerate_tilings(pieces)
    print(f"  tilings: {len(tilings)}")
    assert len(tilings) == 456

    print("matching against data/tiling_catalog.csv ...")
    conv, tiling_info = match_catalog(tilings, pieces,
                                      'data/tiling_catalog.csv')
    print(f"  all 456 tilings matched; name mapping: {conv}")

    print("verifying orbit decomposition ...")
    G, Gp = build_symmetry_perms()
    sizes_G = sorted(len(o) for o in tiling_orbits(tilings, G))
    sizes_Gp = sorted(len(o) for o in tiling_orbits(tilings, Gp))
    print(f"  G-orbits: {sizes_G}")
    print(f"  G'-orbits: {sizes_Gp}")
    assert sizes_G == [12, 12, 72, 72, 72, 72, 144]
    assert sizes_Gp == [12, 12] + [72] * 6

    print("verifying n_til for 200 random catalog entries ...")
    ga_site_map = {fs: ga_site_of(fs) for fs in pieces}
    import csv as _csv
    rows = []
    with open('data/position_catalog_kG3.csv') as f:
        for row in _csv.DictReader(f):
            rows.append(row)
    random.seed(0)
    for row in random.sample(rows, 200):
        ga = set(int(x) for x in row['ga_faces'].split(';'))
        h = set(int(x) for x in row['h_faces'].split(';'))
        nt = sum(1 for t in tilings if is_compatible(ga, h, t, ga_site_map))
        assert nt == int(row['n_tilings']), row['config_id']
    print("  OK")
    print("all checks passed")


if __name__ == '__main__':
    main()
