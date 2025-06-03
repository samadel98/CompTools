#!/usr/bin/env python3
"""
Compute how many solvent (e.g. toluene) molecules fill
a cubic simulation box of given side length, with command-line flags.
"""

import argparse

# Constants
AVOGADRO = 6.02214076e23        # molecules per mol
DEFAULT_DENSITY = 0.867         # g / mL (toluene)
DEFAULT_MOLAR_MASS = 92.14      # g / mol (toluene)
NM_TO_CM = 1e-7                 # 1 nm = 1e-7 cm
CM3_TO_ML = 1.0                 # 1 cm³ = 1 mL

def compute_box_fill(side_nm: float,
                     density: float,
                     molar_mass: float):
    """
    Given a cubic box side in nm, and a liquid’s density (g/mL)
    & molar mass (g/mol), returns
      - box_vol_mL: volume of the box in mL
      - n_moles:      moles of liquid that fill that volume
      - n_molecules:  corresponding number of molecules
    """
    side_cm = side_nm * NM_TO_CM
    box_vol_cm3 = side_cm**3
    box_vol_mL = box_vol_cm3 * CM3_TO_ML

    mass_g      = density * box_vol_mL    # mass in grams
    n_moles     = mass_g / molar_mass     # mol
    n_molecules = n_moles * AVOGADRO

    return box_vol_mL, n_moles, n_molecules

def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate how many solvent molecules fill a cubic box."
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
        '-s', '--box_size',
        type=float,
        required=True,
        help="Box side length in nm"
    )
    return parser.parse_args()

def main():
    args = parse_args()

    side_nm    = args.box_size
    density    = args.density
    molar_mass = args.molarmass

    vol_mL, moles, molecules = compute_box_fill(side_nm, density, molar_mass)

    print(f"Simulation box: {side_nm:.2f} nm per side")
    print(f" → Volume       : {vol_mL:.3e} mL")
    print(f"Liquid parameters:")
    print(f"   • Density    : {density:.3f} g/mL")
    print(f"   • Molar mass : {molar_mass:.3f} g/mol")
    print(f" → Moles        : {moles:.3e} mol")
    print(f" → Molecules    : {molecules:.3e} ≃ {int(molecules):,} molecules")

if __name__ == "__main__":
    main()

