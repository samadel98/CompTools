#!/usr/bin/env python3
"""
Compute how many solvent molecules fill a box (cubic, square‐base, or fully orthogonal)
of given side lengths.
"""

import argparse

# Constants
AVOGADRO = 6.02214076e23    # molecules per mol
DEFAULT_DENSITY = 0.867     # g/mL (e.g. toluene)
DEFAULT_MOLAR_MASS = 92.14  # g/mol (e.g. toluene)
NM_TO_CM = 1e-7             # 1 nm = 1e-7 cm
CM3_TO_ML = 1.0             # 1 cm³ = 1 mL

def compute_box_fill(box_dims_nm: tuple[float, float, float],
                     density: float,
                     molar_mass: float):
    """
    Given box dimensions (x, y, z) in nm, and a liquid’s density (g/mL)
    & molar mass (g/mol), return:
      - box_vol_mL: volume of the box in mL
      - n_moles:    moles of liquid that fill that volume
      - n_molecules: corresponding number of molecules
    """
    # convert each nm → cm
    dims_cm = [dim_nm * NM_TO_CM for dim_nm in box_dims_nm]
    # volume in cm³
    box_vol_cm3 = dims_cm[0] * dims_cm[1] * dims_cm[2]
    # 1 cm³ = 1 mL
    box_vol_mL = box_vol_cm3 * CM3_TO_ML

    mass_g      = density * box_vol_mL     # mass in grams
    n_moles     = mass_g / molar_mass      # mol
    n_molecules = n_moles * AVOGADRO

    return box_vol_mL, n_moles, n_molecules

def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate how many solvent molecules fill a box."
    )
    parser.add_argument(
        '-d', '--density',
        type=float,
        default=DEFAULT_DENSITY,
        help=f"Liquid density in g/mL"
    )
    parser.add_argument(
        '-mm', '--molarmass',
        type=float,
        default=DEFAULT_MOLAR_MASS,
        help=f"Liquid molar mass in g/mol"
    )
    parser.add_argument(
        '-b', '--box_dims',
        type=float,
        nargs='+',
        required=True,
        help=(
            "Box dimension(s) in nm. Provide:\n"
            "  • 1 value → cubic box (all sides = value)\n"
            "  • 2 values → square‐base box (x = y = first, z = second)\n"
            "  • 3 values → orthogonal box (x, y, z)\n"
        )
    )
    return parser.parse_args()

def main():
    args = parse_args()

    raw_dims = args.box_dims
    if len(raw_dims) == 1:
        # cubic: x = y = z = raw_dims[0]
        dims_nm = (raw_dims[0], raw_dims[0], raw_dims[0])
        box_type = "cubic"
    elif len(raw_dims) == 2:
        # square‐base: x = y = raw_dims[0], z = raw_dims[1]
        dims_nm = (raw_dims[0], raw_dims[0], raw_dims[1])
        box_type = "square‐base orthogonal"
    elif len(raw_dims) == 3:
        # fully orthogonal: x, y, z
        dims_nm = tuple(raw_dims)
        box_type = "orthogonal"
    else:
        raise SystemExit("Error: --box_dims must have 1, 2, or 3 values (in nm).")

    density    = args.density
    molar_mass = args.molarmass

    vol_mL, moles, molecules = compute_box_fill(dims_nm, density, molar_mass)

    # Print header with box type and dimensions
    if box_type == "cubic":
        side = dims_nm[0]
        print(f"Simulation box: {side:.2f} nm per side (cubic)")
    elif box_type == "square‐base orthogonal":
        x = dims_nm[0]
        z = dims_nm[2]
        print(
            f"Simulation box: square base {x:.2f}×{x:.2f} nm, "
            f"height {z:.2f} nm (square‐base orthogonal)"
        )
    else:  # orthogonal
        x, y, z = dims_nm
        print(
            f"Simulation box: {x:.2f} nm × {y:.2f} nm × {z:.2f} nm (orthogonal)"
        )

    print(f" → Volume       : {vol_mL:.3e} mL")
    print(f"Liquid parameters:")
    print(f"   • Density    : {density:.3f} g/mL")
    print(f"   • Molar mass : {molar_mass:.3f} g/mol")
    print(f" → Moles        : {moles:.3e} mol")
    print(f" → Molecules    : {molecules:.3e} ≃ {int(molecules):,} molecules")

if __name__ == "__main__":
    main()
