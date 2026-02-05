# VASP post-processing pipeline

Progressive pipeline for analysing VASP (DFT) results. Extensible for future needs.

## Setup

```bash
pip install -r requirements.txt
```

## OSZICAR energy plot

Plot a selected property (default: **E0**) vs ionic step from VASP `OSZICAR` files.

**Usage:**

```bash
# Plot E0 vs step from folder 0_Potential_Basic (figure saved in that folder)
python plot_oszicar.py 0_Potential_Basic

# Custom output path (relative to current directory)
python plot_oszicar.py 0_Potential_Basic -o figures/energy.png

# Plot all OSZICAR files in the folder (OSZICAR, OSZICAR_1, ...) on the same figure
python plot_oszicar.py 0_Potential_Basic --multi

# Plot another property: F, dE, or mag
python plot_oszicar.py 0_Potential_Basic -p F

# Plot energy per atom (requires POSCAR in the same folder)
python plot_oszicar.py 0_Potential_Basic -u per_atom
```

**Options:**

| Option | Default | Description |
|--------|---------|-------------|
| `folder` | `.` | Folder containing OSZICAR (e.g. `0_Potential_Basic`) |
| `--multi` | off | Plot all OSZICAR, OSZICAR_1, … in the folder |
| `-p`, `--property` | E0 | Property: `E0`, `F`, `dE`, `mag` |
| `-u`, `--units` | total | `total` or `per_atom` (per_atom needs POSCAR in folder) |
| `-o`, `--output` | folder/energy_vs_step.png | Output figure path |

Figures are saved at 300 dpi and are suitable for papers/presentations.
