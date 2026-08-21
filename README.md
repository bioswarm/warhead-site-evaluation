# Warhead Site Rank

PyMOL scanner for covalent-warhead site selection. Given a protein complex, it finds **target-interface nucleophiles** and ranks **binder** residues that can be mutated to Cys for linker/warhead installation.

This is a structural heuristic. It does **not** predict reaction rate, covalent yield, or in vivo effect.

---

## Final version (use this)

| | |
|---|---|
| **Script** | [`warhead_site_rank.py`](warhead_site_rank.py) |
| **Version** | **v3.1.3** |
| **Location** | repository root |

The two earlier scripts are archived and are **not** the default entry point:

| Archive path | Notes |
|--------------|-------|
| [`legacy/v2/pymol_warhead_site_rank_v2_1.py`](legacy/v2/pymol_warhead_site_rank_v2_1.py) | v2 line: CB anchor + truncated SASA |
| [`legacy/v3.0.1/warhead_site_rank_v3_0_1.py`](legacy/v3.0.1/warhead_site_rank_v3_0_1.py) | v3.0.1: chemistry profiles + Cys rotamers |

The final script keeps the v3.0.1 architecture, ports residue aliases / object resolution / readable CSV columns from v2, and adds interface gating for existing Cys, ligand clash/path checks, and short PyMOL labels.

---

## Requirements

- PyMOL 2.x / 3.x (`run` to load the script)
- No extra pip packages

## Usage

```pymol
load complex.cif, complex
run /path/to/warhead_site_rank.py
warhead_help
warhead_profiles
warhead_scan_prompt
warhead_scan_prompt advanced=1
```

Direct scan:

```pymol
warhead_scan object_name=complex, target_chain=B, binder_chain=A, \
    profile=fluorosulfate, interface_cutoff=5, candidate_cutoff=12, \
    ideal_min=5, ideal_max=10, label_top=3, out_prefix=warhead_scan
```

Multi-chain or arbitrary selections:

```pymol
warhead_scan object_name=complex, binder_chain=H+L, target_chain=B
warhead_scan object_name=complex, \
    target_sel="chain B+D and resi 1-200", \
    binder_sel="chain A+C", profile=fluorosulfate
```

Force-field residue names:

```pymol
warhead_scan ... target_types=CYX+HIE
warhead_scan ... resn_aliases=CSO:CYS
```

Batch:

```pymol
warhead_batch manifest=/path/jobs.csv, output_dir=/path/results
```

## Chemistry profiles

| profile | Reactive atoms |
|---------|----------------|
| `fluorosulfate` | TYR:OH, LYS:NZ, HIS:ND1/NE2 |
| `broad_nucleophile` | CYS:SG, LYS:NZ, TYR:OH, HIS:ND1/NE2 |
| `cys_target` | CYS:SG |
| `custom` | `reactive_atoms="CYS:SG;LYS:NZ;HIS:ND1\|NE2"` |

## Two primary distances

- `interface_cutoff`: target reactive atom to any binder heavy atom; defines whether a nucleophile sits at the interface (default 5 Å; **not** a chemical law).
- `candidate_cutoff`: modeled Cys-SG to the target reactive atom; search radius for installation sites (default 12 Å).
- `ideal_min` / `ideal_max`: preferred SG distance for the current linker/warhead family.

Existing Cys on the binder are χ1-sampled only if the residue itself is on the protein–protein interface (heavy-atom distance to the target protein ≤ `interface_cutoff`).

## Output

Each scan writes:

- `<prefix>_targets.csv`
- `<prefix>_pairs.csv`
- `<prefix>_binder_rank.csv`
- `<prefix>_metadata.json`

PyMOL labels at most three binder sites:

- binder: `#1 A/Ser74  8.55Å`
- target: `B/Tyr134`
- class is color only: A green / B cyan / C magenta / D gray

Class meaning: A = inspect first; B = first-round candidate; C = backup; D = generally skip first round.

## Public commands

`warhead_help` · `warhead_profiles` · `warhead_scan` · `warhead_scan_prompt` · `warhead_batch`

## Limitations

The script does not enumerate a full maleimide-linker-warhead ensemble. Top hits still need manual rotamer/linker/path inspection, hotspot review, WT vs Cys-mutant vs conjugated binding, and mass-spec confirmation of the cross-link.
