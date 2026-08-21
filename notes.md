# Warhead Site Rank v3.1.3

Final deliverable script: `warhead_site_rank.py`.

This release keeps the reusable v3 architecture (profiles, modeled Cys
rotamers, class caps, Qt prompt, batch/metadata) and lightly ports the most
useful engineering pieces from the v2 script.

## Public commands

```pymol
run /path/to/warhead_site_rank.py
warhead_help
warhead_profiles
warhead_scan_prompt
warhead_scan_prompt advanced=1
warhead_scan ...
warhead_batch manifest=/path/jobs.csv, output_dir=/path/results
```

PyMOL binary used for validation:

```text
/mnt/local_disk1/software/miniforge3/envs/pymol/bin/pymol
```

## Main capabilities

- Accepts chain IDs or arbitrary PyMOL target/binder selections, including multi-chain selections.
- Built-in `fluorosulfate`, `broad_nucleophile`, and `cys_target` chemistry profiles plus custom reactive-atom definitions.
- Explicit interface / candidate / preferred / contact / polar / salt-bridge / internal-contact / clash / path thresholds.
- Builds canonical Cys χ1 rotamers and evaluates modeled SG distance, orientation, clash clearance, solvent-accessibility proxy, and straight-path obstruction.
- Removes the original side chain before modeled SG exposure analysis.
- Precomputes binder-site interface features once and reuses them across target pairs.
- Separates `raw_score`, `uncapped_class`, `final_class`, `class_cap`, and `downgrade_reasons`.
- Natural residue sorting and coordinate-based peptide connectivity for terminal labels.
- Labels at most three binder sites as `#1 A/Ser74  8.55Å` plus a short target tag `B/Tyr134`. Class is color only. `label_top=0`, target-only hits, and no-hit cases are safe.
- CSV/JSON metadata output and manifest-driven batch analysis.

## v3.1.3

- Binder existing Cys are χ1-sampled only when the residue is on the protein-protein interface: heavy-atom distance to the target chain ≤ `interface_cutoff`.
- Target Cys were already interface-gated by the nucleophile `interface_cutoff`.
- Reach prefilter for Cys uses CB or native SG, so an interface Cys whose SG currently points away is not dropped.

## v3.1.2

- Existing Cys keeps the native SG and also tries the three canonical χ1 rotamers (near-duplicates of native are dropped).
- Clash, straight-path obstruction, and bound SG exposure include organic ligands, metals, and inorganic ions from the same object; crystal water is excluded.
- `best_target` / `target_options` use readable tags such as `B/Tyr134 OH`.

## v3.1.1 ranking / robustness fixes

- N/C termini and fragment ends are class-capped to C so they do not look like first-round sites just because geometry is close.
- Preferred-distance flags now use the modeled SG (scored) distance, not the unmutated anchor.
- `class_cap` is blank unless the class was actually lowered.
- Missing chain IDs print the chains that are actually in the object.
- Contact / polar / salt-bridge counts use a spatial grid with the same cutoffs.

## v3.1.0 additions (from v2 + polish)

- Residue aliases: built-in `HIE`/`HID`/`CYX`/`MSE`/terminal `NALA`/`CGLY` mapping, optional `resn_aliases=`, unknown-residue warnings. Scoring uses canonical names; CSV keeps both `*_resn` and `*_resn_std`.
- Robust object resolution: `get_object_list` None/exception normalization and clearer missing-object errors.
- Readable labels/CSV: `A/Ser74`, `A:S74C`, `binder_site` / `mutation` / `target_site`, richer PyMOL labels.
- Chain-overlap check when both sides are specified as chain IDs.

## Example

```pymol
load complex.cif, complex
run /path/to/warhead_site_rank.py
warhead_scan object_name=complex, target_chain=B, binder_chain=A, \
    profile=fluorosulfate, interface_cutoff=5, candidate_cutoff=12, \
    ideal_min=5, ideal_max=10, label_top=3, out_prefix=warhead_scan
```

## Design / plan docs

- `docs/superpowers/specs/2026-08-21-warhead-site-rank-design.md`
- `docs/superpowers/plans/2026-08-21-warhead-site-rank.md`
