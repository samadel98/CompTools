#!/usr/bin/env python3

"""
delete_workflow.py

This script is designed to help users interactively edit atomistic structures stored in .xyz files, 
commonly used in computational chemistry and materials science. It allows you to:

    1. Replace the nearest neighbor atoms of a given type (e.g., substitute the closest 'Cs' near 100 randomly 
       selected 'Br' atoms by 'Rb').
    2. Remove the nearest neighbor atom of a given type for randomly selected reference atoms, and then 
       convert those reference atoms to another type (e.g., delete the closest 'Br' neighbor to 100 randomly 
       selected 'Br' atoms, and change those 100 atoms to 'F').
       
You can choose to perform only the first step, or both steps, and will be prompted after each main step.
The script prints a summary of changes and asks for confirmation before making modifications.

Usage:
    python delete_workflow.py <input_xyz_file> <output_xyz_file>

Example:
    python delete_workflow.py input.xyz output.xyz

Author: (Your Name)
Date: (Date)
"""

import sys
import random
import math

# --- Utility Functions (Unchanged) ---

def read_xyz(filename):
    """
    Reads an .xyz file and returns:
      - num_atoms (int)
      - comment (str)
      - atoms_list (list of tuples): (atom_type, x, y, z)
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    num_atoms = int(lines[0].strip())
    comment = lines[1].strip()
    
    atoms_list = []
    for i in range(2, 2 + num_atoms):
        parts = lines[i].split()
        atom_type = parts[0]
        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        atoms_list.append((atom_type, x, y, z))
    
    return num_atoms, comment, atoms_list


def write_xyz(filename, atoms_list, comment=""):
    """
    Writes an .xyz file given a list of atoms in the format:
      [(atom_type, x, y, z), ...]
    """
    with open(filename, 'w') as f:
        f.write(f"{len(atoms_list)}\n")
        f.write(comment + "\n")
        for atom in atoms_list:
            atom_type, x, y, z = atom
            f.write(f"{atom_type} {x:.6f} {y:.6f} {z:.6f}\n")


def distance(atom_a, atom_b):
    """
    Returns the Euclidean distance between two atoms:
    Each atom is (atom_type, x, y, z).
    """
    _, ax, ay, az = atom_a
    _, bx, by, bz = atom_b
    return math.sqrt((ax - bx)**2 + (ay - by)**2 + (az - bz)**2)


def find_closest_atom(
    reference_atom,
    atoms_list,
    target_type,
    already_marked=None,
    ref_idx=None,
    max_attempts=10
):
    """
    Finds the single closest atom in `atoms_list` of type `target_type`,
    skipping:
      - any atom whose index is in 'already_marked' (if provided),
      - the reference index (ref_idx),
      - any atom that is None in the list.
    
    Returns the index of the closest valid atom or None if none found.
    """
    if already_marked is None:
        already_marked = set()
    
    candidates = []
    for i, atom in enumerate(atoms_list):
        if atom is None:
            continue  # skip deleted/None
        if ref_idx is not None and i == ref_idx:
            continue  # don't let the reference see itself as neighbor
        if i in already_marked:
            continue
        
        # Must match target type
        if atom[0] == target_type:
            dist = distance(reference_atom, atom)
            candidates.append((i, dist))
    
    # Sort by ascending distance
    candidates.sort(key=lambda x: x[1])
    
    for attempt_i, (idx, dist) in enumerate(candidates):
        if attempt_i >= max_attempts:
            break
        if idx not in already_marked:
            return idx
    
    return None


def count_atom_types(atoms_list):
    """
    Returns a dictionary of atom_type -> count in the given atoms_list (skips None).
    """
    type_counts = {}
    for atom in atoms_list:
        if atom is None:
            continue
        a_type, x, y, z = atom
        type_counts[a_type] = type_counts.get(a_type, 0) + 1
    return type_counts


def confirm_step(message):
    """
    Asks user (y/n) to confirm if they'd like to proceed.
    If 'n', then exits. If 'y', continues.
    """
    while True:
        choice = input(f"{message} (y/n): ").strip().lower()
        if choice == 'y':
            return True
        elif choice == 'n':
            print("Operation cancelled by user. Exiting.")
            sys.exit(0)
        else:
            print("Please type 'y' for yes or 'n' for no.")

def proceed_or_save(atoms_list, comment, output_xyz, step_descr):
    """
    After each main step, ask the user if they want to proceed to the next step or save and exit.
    """
    while True:
        choice = input(
            f"\nWould you like to proceed to the next step, or save and exit now?\n"
            f"Type 'c' to continue, 's' to save and exit: ").strip().lower()
        if choice == 'c':
            return True
        elif choice == 's':
            write_xyz(output_xyz, [atom for atom in atoms_list if atom is not None], comment)
            print(f"\n{step_descr} complete. Changes saved to '{output_xyz}'. Exiting now.\n")
            sys.exit(0)
        else:
            print("Please type 'c' to continue or 's' to save and exit.")


# --- Main Workflow ---

def main():
    if len(sys.argv) != 3:
        print(
            f"Usage: python {sys.argv[0]} <input_xyz_file> <output_xyz_file>\n"
            "For example: python delete_workflow.py input.xyz output.xyz"
        )
        sys.exit(1)
    
    input_xyz = sys.argv[1]
    output_xyz = sys.argv[2]
    
    print("\n===== Atomistic Structure Modification Workflow =====")
    print("This tool allows you to edit .xyz files by replacing or deleting atoms based on proximity.")
    print("You can perform up to two sequential steps, with the option to save and exit after each step.\n")
    
    # Step 0: Read input
    num_atoms, comment, atoms_list = read_xyz(input_xyz)
    initial_counts = count_atom_types(atoms_list)
    print("Initial atom-type counts:")
    for t, c in sorted(initial_counts.items()):
        print(f"  {t} : {c}")
    print()
    
    # --------- Step 1 ---------
    print("--- Step 1: Replace the Closest Neighbor's Type ---")
    print("In this step, you will:")
    print("  - Randomly select a number of reference atoms of a specified type ('type1').")
    print("  - For each, find the closest neighbor atom of a different specified type ('type2').")
    print("  - Change each found neighbor atom to a new type you choose.\n")
    
    type1_for_replacement = input("Enter the symbol of the REFERENCE atom type to pick (e.g., 'Cl'): ").strip()
    n_pick_for_replacement = int(input(f"How many '{type1_for_replacement}' atoms do you want to select for replacement? "))
    type2_to_replace = input("Enter the symbol of the NEIGHBOR atom type you want to replace (e.g., 'Cs'): ").strip()
    type2_new_type = input("Enter the symbol for the NEW atom type to use as replacement (e.g., 'Rb'): ").strip()
    
    current_counts = count_atom_types(atoms_list)
    available_count = current_counts.get(type1_for_replacement, 0)
    print(f"\nFound {available_count} '{type1_for_replacement}' atoms in the structure.")
    print(f"You have requested to randomly select {n_pick_for_replacement} for this operation.")
    
    confirm_step(f"Do you want to proceed with Step 1 (replace neighbors of {type1_for_replacement})?")
    
    # Step 1: Do the replacement
    type1_indices = [i for i, atm in enumerate(atoms_list) if atm is not None and atm[0] == type1_for_replacement]
    if len(type1_indices) < n_pick_for_replacement:
        print(f"Error: Not enough '{type1_for_replacement}' atoms to pick {n_pick_for_replacement}. Only found {len(type1_indices)}.")
        sys.exit(1)
    
    chosen_type1_indices = random.sample(type1_indices, n_pick_for_replacement)
    replaced_indices = set()
    
    for ref_idx in chosen_type1_indices:
        ref_atom = atoms_list[ref_idx]
        if ref_atom is None:
            continue  # Should not happen
        closest_idx = find_closest_atom(
            reference_atom=ref_atom,
            atoms_list=atoms_list,
            target_type=type2_to_replace,
            already_marked=replaced_indices,
            ref_idx=ref_idx,
            max_attempts=10
        )
        if closest_idx is not None:
            old_type, x, y, z = atoms_list[closest_idx]
            atoms_list[closest_idx] = (type2_new_type, x, y, z)
            replaced_indices.add(closest_idx)
    
    after_step1_counts = count_atom_types(atoms_list)
    print("\nAtom counts after Step 1:")
    for t, c in sorted(after_step1_counts.items()):
        print(f"  {t} : {c}")
    
    confirm_step("Does the atom count above match your expectation for Step 1?")
    
    # Prompt: Proceed to next step or save and exit?
    proceed_or_save(atoms_list, comment, output_xyz, "Step 1")
    
    # --------- Step 2 ---------
    print("\n--- Step 2: Delete Closest Neighbor and Convert Reference Atoms ---")
    print("In this step, you will:")
    print("  - Randomly select reference atoms ('type1').")
    print("  - For each, find and DELETE their closest neighbor atom of a specified type.")
    print("  - Then, change the selected reference atom to a new atom type.\n")
    
    type1_for_deletion = input("Enter the symbol of the REFERENCE atom type to pick (e.g., 'Cl'): ").strip()
    n_pick_for_deletion = int(input(f"How many '{type1_for_deletion}' atoms do you want to process in this step? "))
    type_for_neighbor_deletion = input("Enter the symbol of the NEIGHBOR atom type you want to delete (e.g., 'Cl'): ").strip()
    new_type_for_chosen_del = input("Enter the symbol for the NEW atom type to convert the reference atoms to (e.g., 'F'): ").strip()
    
    current_counts2 = count_atom_types(atoms_list)
    available_count2 = current_counts2.get(type1_for_deletion, 0)
    print(f"\nFound {available_count2} '{type1_for_deletion}' atoms in the structure.")
    print(f"You have requested to randomly select {n_pick_for_deletion} atoms for this operation.")
    
    confirm_step(f"Do you want to proceed with Step 2 (delete neighbors of {type1_for_deletion} and convert their type)?")
    
    # Step 2: One by one delete & convert
    for iteration in range(1, n_pick_for_deletion + 1):
        current_type1_indices = [i for i, atm in enumerate(atoms_list)
                                 if atm is not None and atm[0] == type1_for_deletion]
        if not current_type1_indices:
            print("No more atoms of the selected type remain to pick from. Stopping early.")
            break
        
        chosen_ref_idx = random.choice(current_type1_indices)
        ref_atom = atoms_list[chosen_ref_idx]  # Still type1_for_deletion
        closest_idx = find_closest_atom(
            reference_atom=ref_atom,
            atoms_list=atoms_list,
            target_type=type_for_neighbor_deletion,
            already_marked=None,
            ref_idx=chosen_ref_idx,
            max_attempts=50
        )
        if closest_idx is not None:
            atoms_list[closest_idx] = None
        
        # Now convert the chosen reference atom to the new type
        old_type, x, y, z = ref_atom
        atoms_list[chosen_ref_idx] = (new_type_for_chosen_del, x, y, z)
    
    final_atoms_list = [atom for atom in atoms_list if atom is not None]
    
    after_step2_counts = count_atom_types(final_atoms_list)
    print("\nAtom counts after Step 2:")
    for t, c in sorted(after_step2_counts.items()):
        print(f"  {t} : {c}")
    
    confirm_step("Does the atom count above match your expectation for Step 2?")
    
    # Last prompt: Save and exit
    write_xyz(output_xyz, final_atoms_list, comment="Modified structure")
    print(f"\nAll steps complete! Final structure written to '{output_xyz}'.\n")

if __name__ == "__main__":
    main()

