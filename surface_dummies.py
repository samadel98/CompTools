#!/usr/bin/env python3
"""
surface_replace_tool.py

A robust, interactive CLI tool to replace specific atom types in a molecular structure (XYZ file)
using SCM PLAMS + CAT's `replace_surface`, with intelligent suggestions for replacement ("dummy") atoms.

Features:
- Shows atom types and counts at every step.
- Lets the user replace atoms by either a **fraction** or an **exact number**.
- For exact number, does a dry run to ensure the right number of replacements, regardless of randomness in the replacement logic.
- Suggests dummy atoms from the same periodic column (using mendeleev, with symbols and names in a table).
- Saves the final molecule to a new XYZ file and displays a summary.

Requirements:
- scm.plams
- CAT.recipes.replace_surface
- mendeleev

USAGE:
    python surface_replace_tool.py <input.xyz> <output.xyz>

Example:
    python surface_replace_tool.py structure.xyz structure_modified.xyz

"""

import sys
import copy
from collections import Counter
from scm.plams import Molecule
from CAT.recipes import replace_surface
from mendeleev import element, get_all_elements

def print_atom_counts(mol, msg=None):
    """Print a formatted summary of atom types and their counts in the molecule."""
    counts = Counter(a.symbol for a in mol)
    if msg:
        print(msg)
    print("\nAtom types and counts in the loaded structure:")
    for symbol, count in sorted(counts.items()):
        print(f"  {symbol:<3}: {count}")
    print()

def same_column_candidates(symbol):
    """
    Returns a list of mendeleev.Element objects from the same group (column) as the input symbol,
    excluding the symbol itself.
    """
    try:
        el = element(symbol)
        if el.group_id is None:
            return []
        group = el.group_id
        # Exclude the symbol itself and long/rare synthetic ones
        candidates = [e for e in get_all_elements() if e.group_id == group and e.symbol != symbol]
        candidates = [e for e in candidates if len(e.symbol) <= 2]
        return candidates
    except Exception:
        return []

def show_candidates_table(candidates):
    """Pretty-prints a table of candidate dummy atom symbols and names."""
    print("\nSuggested dummy atoms from the same periodic column:")
    print(f"{'Symbol':^8} | {'Name':^15}")
    print("-" * 25)
    for el in candidates:
        print(f"{el.symbol:^8} | {el.name:^15}")
    print()

def ask_fraction_or_number(symbol, original_count, mol, symbol_new):
    """
    Asks user if they want to replace by fraction or by number.
    If by number, does a dry run with f=1.0 to calculate the right fraction for exact replacement.
    Returns (f, n_target)
    """
    while True:
        mode = input(f"\nReplace '{symbol}' atoms by (f)raction or (n)umber? [f/n]: ").lower().strip()
        if mode in ('f', 'n'):
            break
        print("Invalid input. Enter 'f' or 'n'.")
    if mode == 'f':
        while True:
            try:
                f = float(input("Enter the fraction f (between 0.0 and 1.0): ").strip())
                if 0 <= f <= 1:
                    break
                else:
                    print("Fraction must be between 0.0 and 1.0.")
            except Exception:
                print("Invalid input, try again.")
        n_target = int(round(f * original_count))
        return f, n_target
    else:  # by number
        while True:
            try:
                n_target = int(input(f"How many '{symbol}' atoms to replace (max {original_count})? ").strip())
                if 0 < n_target <= original_count:
                    break
                else:
                    print(f"Enter a value between 1 and {original_count}.")
            except Exception:
                print("Invalid input, try again.")
        # Dry run: full replacement, count how many can actually be replaced!
        temp_mol = copy.deepcopy(mol)
        temp_mol = replace_surface(
            temp_mol,
            symbol=symbol,
            symbol_new=symbol_new,
            f=1.0,
            mode='uniform',
            displacement_factor=0.7
        )
        dummy_count = sum(1 for a in temp_mol if a.symbol == symbol_new)
        if dummy_count == 0:
            print(f"Error: No '{symbol_new}' atoms created in full-replacement dry run. Aborting.")
            sys.exit(2)
        f_real = n_target / dummy_count
        print(f"  (Internal dry run: Replacing all '{symbol}' → '{symbol_new}' gives {dummy_count} dummy atoms. " +
              f"Applying f = {f_real:.5f} for your request of {n_target})")
        return f_real, n_target

def main():
    # --- Help / usage message ---
    if len(sys.argv) != 3:
        print(__doc__)
        print("ERROR: Must provide <input.xyz> <output.xyz>")
        sys.exit(1)

    old_file, new_file = sys.argv[1:3]
    # --- Load molecule and print initial atom stats ---
    mol = Molecule(old_file)
    print_atom_counts(mol)

    print("Enter your surface replacement instructions.")
    print("Leave 'symbol' empty (just press ENTER) when you are done.\n")

    # Save the original atom counts for correct handling of replacements by number
    original_counts = Counter(a.symbol for a in mol)

    # --- Main interactive loop for atom replacement ---
    while True:
        symbol = input("Enter the atom type to replace [or ENTER to finish]: ").strip()
        if not symbol:
            break

        if symbol not in original_counts or original_counts[symbol] == 0:
            print(f"No atoms of type '{symbol}' found. Please try again.")
            continue

        # Candidates: same periodic column, pretty-printed table
        candidates = same_column_candidates(symbol)
        if candidates:
            show_candidates_table(candidates)
        else:
            print("(No candidates from the same periodic column found.)")

        symbol_new = input(f"Enter the dummy atom type for '{symbol}': ").strip()
        if not symbol_new:
            print("No dummy atom type entered, skipping replacement.")
            continue

        # Choose fraction or number (dry run if needed for exact number)
        f, n_target = ask_fraction_or_number(symbol, original_counts[symbol], mol, symbol_new)
        print(f"\nReplacing {n_target} of '{symbol}' atoms with '{symbol_new}'...")

        # --- Atom replacement ---
        mol = replace_surface(
            mol,
            symbol=symbol,
            symbol_new=symbol_new,
            f=f,
            mode='uniform',
            displacement_factor=0.7
        )
        print("Done.")
        print_atom_counts(mol)

    # --- Save and print final file stats ---
    mol.write(new_file)
    print(f"Final model written to '{new_file}'.\n")

    final_mol = Molecule(new_file)
    print_atom_counts(final_mol, msg="Final atom types and counts in the output file:")

if __name__ == "__main__":
    main()

