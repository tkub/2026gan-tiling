# Tiling decomposition multiplicity of the EC constraint on GaN(0001)

Data accompanying the paper:

> T. Kuboyama, A. Kusaba, K. Kawka, and P. Kempisty, *Tiling decomposition multiplicity
> predicts stability of GaN(0001) surface reconstructions*,
> [arXiv:2607.11105](https://arxiv.org/abs/2607.11105) (2026).
> <!-- Replace with the journal reference upon publication. -->

The paper enumerates all 416,683 EC-compatible adatom configurations
(3 Ga adatoms + 18 H atoms) on the GaN(0001)-(6×6) surface and shows
that the number of rhombus tilings compatible with a configuration
(the tiling decomposition multiplicity, `n_til`) predicts the most
stable configuration within each of the 14 symmetry classes of Ga
placements.

## Conventions

The (6×6) toroidal triangular lattice has 72 triangular faces indexed
by `face_id = 12q + 2r + t` with axial coordinates `q, r ∈ {0,…,5}`
and sublattice index `t`:

- `t = 0`: downward-pointing faces (N sublattice). Ga adatom sites.
- `t = 1`: upward-pointing faces (Ga sublattice). H adsorption sites.

Face lists in the CSV files are semicolon-separated `face_id` values.
Symmetry classification uses the sublattice-preserving group
G′ = T⋊D₃ (|G′| = 216). Ga placement classes are numbered 1–14 in the
paper (CSV column `ga_orbit` is 0-indexed: class = `ga_orbit` + 1).
Tiling classes are labeled T1–T7a/T7b in the paper; the column
`class_id` in `tiling_catalog.csv` gives the classification under the
full group G (T7a/T7b merged as class 7). The legacy orientation
names `up`, `ne`, `nw` in `tiling_catalog.csv` correspond to
orientations 1, 2, 3 of the paper.

## Data files (`data/`)

| File | Rows | Description |
|:---|---:|:---|
| `position_catalog_kG3.csv` | 416,683 | Complete catalog of EC-compatible configurations. Columns: `config_id`, `ga_orbit` (0-indexed class), `ga_faces`, `h_faces`, `n_tilings` (= n_til). |
| `tiling_catalog.csv` | 456 | All rhombus tilings of the bare lattice. Columns: `global_id`, `class_id` (G-classification), `is_visual_rep`, `r1`–`r9` (piece orientation and anchor), `n_up`, `n_ne`, `n_nw` (orientation counts n₁, n₂, n₃). |
| `ga_orbit_summary_kG3.csv` | 14 | The 14 Ga placement classes: `orbit_id` (0-indexed), `orbit_size` (\|O\|), `canonical_ga`, `d_sq_signature` (squared pair distances of the Ga triangle), `n_patterns`. |
| `relax_results.csv` | 2,529 | MLIP (UMA) relaxation energies of the stratified sample. Columns: `config_id`, `n_tilings`, `energy_eV` (total energy), `converged`. |
| `dft_results.csv` | 25 | Cross-class DFT validation set (all 24 n_til-max candidate configurations plus the Class-2 minimum, spanning the 14 Ga placement classes): `config_id`, `E_dft_total_eV`. |
| `dft_uma_mapped.csv` | 685 | Fixed-placement DFT validation set (Class-7 Ga placement): `dataset` (sampling run of the source database, `energy-1st`/`energy-2nd`), `structure_id` (id in the source database), `config_id` (−1 for EC-incompatible structures), `n_til` (0 for EC-incompatible), `dE_DFT_eV`, `dE_UMA_eV` (relative energies), `ga_faces`, `h_faces` (adatom positions in the face-index convention of `position_catalog_kG3.csv`; used, e.g., to compute the local frustration count `n_adj`). The structures and DFT energies originate from the two Bayesian-optimization runs of K. Kawka *et al.*, J. Appl. Phys. **135**, 225302 (2024). |
| `relax_needed_nadj0.csv` | 76 | The n_adj = 0 configurations outside the stratified sample (config_id, ga_orbit, n_tilings, ga_faces, h_faces); together with the 26 members of `relax_results.csv` they form the complete n_adj = 0 set (102). |
| `relax_nadj0_untileable.csv` | 88 | Bare-site arrangements with n_adj = 0 that admit no compatible tiling, at the 14 canonical Ga placements (pseudo-ids 900000+; same columns). |
| `relax_all167_out.csv` | 167 | MLIP relaxation results (identical protocol to `relax_results.csv`) for 3 calibration configurations, the 76 configurations above, and the 88 untileable arrangements. |

Energies are per (6×6) surface cell. MLIP energies were computed with
the UMA machine-learning interatomic potential via ASE; DFT settings
are described in the paper (Sec. II E).

## Scripts (`scripts/`)

Python ≥ 3.8, standard library only. Run from the repository root,
in this order:

```
python3 scripts/enumeration.py
python3 scripts/compute_descriptors_relaxed.py
python3 scripts/paper_statistics.py
```

| Script | Purpose |
|:---|:---|
| `enumeration.py` | Library for the tiling enumeration and symmetry analysis. Run as a script, it re-derives everything from the definitions and verifies the shipped catalogs: the 456 tilings (exact set cover), the orbit decompositions under G (7 classes) and G′ (8 classes), and sampled `n_til` values. |
| `compute_descriptors_relaxed.py` | Recomputes `n_til` (verified against the catalog), the geometric descriptors of the paper (I_iso and its Shannon variant, V_HH, mean Ga–H distance), and the local frustration count `n_adj` (number of nearest-neighbor pairs of bare sites) for the 2,529 relaxed configurations, and writes `data/descriptors_relaxed.csv`. |
| `paper_statistics.py` | Prints the statistics quoted in the paper: catalog statistics, the per-class verification of the n_til-max rule (Table III), descriptor correlations, and the MLIP–DFT comparison for the two DFT validation sets. |
| `verify_nadj0_arrangements.py` | Verifies the bare-site arrangement counts of Sec. III F: 246 arrangements with n_adj = 0 at the 14 canonical placements, of which 88 are untileable, and the reduction of the 158 tileable ones to the 102-configuration census. |
| `check_ising_model.py` | Re-evaluates the data-driven Ising model of Kawka et al. (2024) over the 14,896 configurations at the Class-7 placement and confirms its unique minimum, configuration 109117. |

## License

The code (`scripts/`) is released under the MIT License (see
`LICENSE`). The data (`data/`) are released under the Creative
Commons Attribution 4.0 International License (CC BY 4.0; see
`LICENSE-DATA`). If you use these data, please cite the paper above
(see `CITATION.cff`).

## Contact

Tetsuji Kuboyama (Gakushuin University), kuboyama@tk.cc.gakushuin.ac.jp
