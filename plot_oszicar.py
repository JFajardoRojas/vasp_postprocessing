#!/usr/bin/env python3
"""
VASP OSZICAR post-processing: plot energy (or other property) vs ionic step.

Reads OSZICAR files from a given folder and produces a research-quality plot
of the selected property (default: E0) as a function of the ionic cycle (step).
Supports multiple OSZICAR files (OSZICAR, OSZICAR_1, ...) when --multi is set.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless-friendly: save to file, no display required
import matplotlib.pyplot as plt
import numpy as np

# Summary line format: "   N F= value E0= value  d E = value  mag= value"
# Example: "   1 F= -.13179597E+03 E0= -.13180749E+03  d E =-.131796E+03  mag=    -0.0279"
# Numbers can be scientific (E+03, E-02) or plain floats; allow optional minus in exponent.
SUMMARY_LINE = re.compile(
    r"^\s*(\d+)\s+F=\s*([-\d.E+-]+)\s+E0=\s*([-\d.E+-]+)\s+d E =\s*([-\d.E+-]+)\s+mag=\s*([-\d.E+-]+)\s*$"
)

PROPERTY_NAMES = ("E0", "F", "dE", "mag")
PROPERTY_INDEX = {"E0": 2, "F": 1, "dE": 3, "mag": 4}
PROPERTY_LABELS = {
    "E0": r"$E_0$ (eV)",
    "F": r"$F$ (eV)",
    "dE": r"$\Delta E$ (eV)",
    "mag": "mag",
}

# For per-atom labels
PROPERTY_LABELS_PER_ATOM = {
    "E0": r"$E_0$ (eV/atom)",
    "F": r"$F$ (eV/atom)",
    "dE": r"$\Delta E$ (eV/atom)",
    "mag": "mag/atom",
}


def count_atoms_poscar(folder: Path) -> int:
    """Return number of atoms from POSCAR: count coordinate rows after 'direct'/'Cartesian'."""
    poscar = folder / "POSCAR"
    if not poscar.is_file():
        raise FileNotFoundError(f"POSCAR not found in {folder} (required for --units per_atom)")
    # Match a line that looks like "direct" or "Cartesian" (case-insensitive)
    coord_keyword = re.compile(r"^\s*(direct|cartesian)\s*$", re.I)
    # Line with at least 3 floats (x y z coordinates)
    coord_line = re.compile(r"^\s*[-\d.E+]+[\s\d.E+-]+[-\d.E+]")
    n_atoms = 0
    after_keyword = False
    with open(poscar) as f:
        for line in f:
            if coord_keyword.match(line.strip()):
                after_keyword = True
                continue
            if after_keyword:
                stripped = line.strip()
                if not stripped:
                    break  # blank line ends coordinate block
                parts = stripped.split()
                if len(parts) >= 3:
                    try:
                        float(parts[0])
                        float(parts[1])
                        float(parts[2])
                        n_atoms += 1
                    except ValueError:
                        break
                else:
                    break
    if n_atoms == 0:
        raise ValueError(f"Could not find coordinate block or atoms in {poscar}")
    return n_atoms


def parse_oszicar(path: Path, property_name: str) -> tuple[np.ndarray, np.ndarray]:
    """Parse an OSZICAR file and return (cycles, values) for the given property."""
    cycles = []
    values = []
    idx = PROPERTY_INDEX[property_name]
    with open(path) as f:
        for line in f:
            m = SUMMARY_LINE.match(line)
            if m:
                cycles.append(int(m.group(1)))
                values.append(float(m.group(idx)))
    return np.array(cycles), np.array(values)


def find_oszicar_files(folder: Path, multi: bool) -> list[Path]:
    """Return OSZICAR file paths to use. If multi, all OSZICAR*; else only OSZICAR."""
    folder = folder.resolve()
    if not folder.is_dir():
        raise FileNotFoundError(f"Folder not found: {folder}")

    # Collect OSZICAR and OSZICAR_N (N = 1, 2, ...)
    candidates: list[tuple[int, Path]] = []
    for p in folder.iterdir():
        if p.name == "OSZICAR" and p.is_file():
            candidates.append((0, p))
        elif p.name.startswith("OSZICAR_") and p.is_file():
            try:
                n = int(p.name.split("_", 1)[1])
                candidates.append((n, p))
            except ValueError:
                pass
    candidates.sort(key=lambda x: x[0])

    if not candidates:
        raise FileNotFoundError(f"No OSZICAR file found in {folder}")

    if multi:
        return [p for _, p in candidates]
    return [candidates[0][1]]


def plot(
    data: list[tuple[str, np.ndarray, np.ndarray]],
    property_name: str,
    out_path: Path,
    per_atom: bool = False,
) -> None:
    """Build research-quality matplotlib figure."""
    plt.rcParams["font.size"] = 11
    plt.rcParams["axes.labelsize"] = 12
    plt.rcParams["axes.titlesize"] = 13
    plt.rcParams["legend.fontsize"] = 10
    plt.rcParams["figure.dpi"] = 150

    fig, ax = plt.subplots(figsize=(6, 4))
    labels = PROPERTY_LABELS_PER_ATOM if per_atom else PROPERTY_LABELS
    ylabel = labels.get(property_name, property_name)

    for label, cycles, values in data:
        ax.plot(cycles, values, "-o", markersize=3, linewidth=1.2, label=label)

    ax.set_xlabel("Ionic step")
    ax.set_ylabel(ylabel)
    ax.grid(True, linestyle="--", alpha=0.7)
    ax.legend(loc="best")
    fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"Figure saved: {out_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot VASP OSZICAR property vs ionic step (research-quality).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "folder",
        type=Path,
        nargs="?",
        default=Path("."),
        help="Folder containing OSZICAR (e.g. 0_Potential_Basic)",
    )
    parser.add_argument(
        "--multi",
        action="store_true",
        help="Plot all OSZICAR files (OSZICAR, OSZICAR_1, ...) when present",
    )
    parser.add_argument(
        "-p",
        "--property",
        choices=PROPERTY_NAMES,
        default="E0",
        metavar="PROP",
        help="Property to plot: E0, F, dE, mag",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output figure path (e.g. energy_vs_step.png)",
    )
    parser.add_argument(
        "-u",
        "--units",
        choices=("total", "per_atom"),
        default="total",
        metavar="UNITS",
        help="Units: total (default) or per_atom (requires POSCAR in folder)",
    )
    args = parser.parse_args()

    n_atoms: int | None = None
    if args.units == "per_atom":
        n_atoms = count_atoms_poscar(args.folder.resolve())

    paths = find_oszicar_files(args.folder, args.multi)
    data = []
    for p in paths:
        label = p.name
        cycles, values = parse_oszicar(p, args.property)
        if cycles.size == 0:
            print(f"Warning: no summary lines in {p}")
            continue
        if n_atoms is not None:
            values = values / n_atoms
        data.append((label, cycles, values))

    if not data:
        raise SystemExit("No data to plot.")

    out = args.output
    if out is None:
        out = args.folder.resolve() / "energy_vs_step.png"
    else:
        out = Path(out)
        if not out.is_absolute():
            out = Path.cwd() / out
    plot(data, args.property, out, per_atom=(args.units == "per_atom"))


if __name__ == "__main__":
    main()
