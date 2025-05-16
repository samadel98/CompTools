import re
import os
import matplotlib.pyplot as plt

# ------------------ Utility Functions ------------------

def find_minimum_error_line(log_file):
    """
    Reads the log file and searches for lines containing 'error = <value>'.
    Returns a tuple (min_error_value, line_index, all_lines) where:
      - min_error_value: smallest error (float)
      - line_index: index of that line
      - all_lines: list of all lines in the log file.
    """
    min_error = float('inf')
    min_line_index = None
    min_error_value = None
    lines = []
    
    with open(log_file, 'r') as f:
        lines = f.readlines()
    
    for i, line in enumerate(lines):
        match = re.search(r"error\s*=\s*([-\d\.Ee]+)", line)
        if match:
            error_val = float(match.group(1))
            if error_val < min_error:
                min_error = error_val
                min_line_index = i
                min_error_value = error_val
                
    return min_error_value, min_line_index, lines

def find_folder_for_line(lines, line_index):
    """
    Given the log lines and a line index (where the error was found),
    searches upward for the last occurrence of a line with "JOB md_settings.<number>".
    Returns the folder number as a string (preserving any leading zeros).
    """
    for i in range(line_index, -1, -1):
        match = re.search(r"JOB\s+md_settings\.(\d+)", lines[i])
        if match:
            return match.group(1)  # e.g., "003"
    return None

def get_best_folder(log_file="armc.log"):
    """
    Uses the above functions to determine the folder corresponding to the minimum error move.
    Returns (folder_str, min_error) where folder_str is e.g. "003".
    """
    min_error, error_line_index, lines = find_minimum_error_line(log_file)
    if error_line_index is None:
        return None, None
    folder_str = find_folder_for_line(lines, error_line_index)
    return folder_str, min_error

# ------------------ CP2K Input Parsing Functions ------------------

def parse_cp2k_in(cp2k_in_path):
    """
    Parses a CP2K input file for:
      - Atomic charges (in &CHARGE blocks)
      - Lennard-Jones parameters (in &LENNARD-JONES blocks)
    Returns a dictionary with keys 'charges' and 'lj'.
    """
    data = {
        'charges': {},
        'lj': []
    }
    
    in_forcefield = False
    in_charge = False
    in_lj = False
    current_atom = None
    current_pair = []
    
    with open(cp2k_in_path, 'r') as f:
        for line in f:
            stripped = line.strip()
            
            if stripped.upper().startswith("&FORCEFIELD"):
                in_forcefield = True
                continue
            elif stripped.upper().startswith("&END"):
                if in_charge:
                    in_charge = False
                elif in_lj:
                    in_lj = False
                    current_pair = []
                else:
                    in_forcefield = False
                continue
            
            if not in_forcefield:
                continue
            
            if stripped.upper().startswith("&CHARGE"):
                in_charge = True
                continue
            elif stripped.upper().startswith("&LENNARD-JONES"):
                in_lj = True
                continue
            
            if in_charge:
                m_atom = re.search(r"ATOM\s+(\S+)", stripped, re.IGNORECASE)
                if m_atom:
                    current_atom = m_atom.group(1)
                m_charge = re.search(r"CHARGE\s+([-\d\.Ee]+)", stripped, re.IGNORECASE)
                if m_charge and current_atom is not None:
                    data['charges'][current_atom] = float(m_charge.group(1))
            
            if in_lj:
                m_atoms = re.search(r"ATOMS\s+(\S+)\s+(\S+)", stripped, re.IGNORECASE)
                if m_atoms:
                    current_pair = [m_atoms.group(1), m_atoms.group(2)]
                m_epsilon = re.search(r"EPSILON\s+(?:\[.*?\])?\s*([-\d\.Ee]+)", stripped, re.IGNORECASE)
                if m_epsilon and current_pair:
                    epsilon_val = float(m_epsilon.group(1))
                    data['lj'].append({
                        'pair': tuple(current_pair),
                        'epsilon': epsilon_val,
                        'sigma': None
                    })
                m_sigma = re.search(r"SIGMA\s+(?:\[.*?\])?\s*([-\d\.Ee]+)", stripped, re.IGNORECASE)
                if m_sigma and current_pair:
                    sigma_val = float(m_sigma.group(1))
                    if data['lj'] and data['lj'][-1]['pair'] == tuple(current_pair) and data['lj'][-1]['sigma'] is None:
                        data['lj'][-1]['sigma'] = sigma_val
                    else:
                        data['lj'].append({
                            'pair': tuple(current_pair),
                            'epsilon': None,
                            'sigma': sigma_val
                        })
    return data

def print_parameters_table(parsed_data):
    """
    Prints atomic charges and LJ parameters in formatted tables.
    Atomic charges are printed with full precision.
    LJ parameters are printed with 7 significant digits.
    """
    print("Atomic Charges:")
    print("-" * 31)
    print("| {:<4} | {:>15} |".format("Atom", "Charge (e)"))
    print("|" + "-"*4 + "+" + "-"*17 + "|")
    for atom, charge in parsed_data['charges'].items():
        # 15 significant digits for charges
        print("| {:<4} | {:>15.15g} |".format(atom, charge))
    print("-" * 31)
    print()
    
    print("Lennard-Jones Parameters:")
    print("-" * 54)
    print("| {:<9} | {:>15} | {:>15} |".format("Atom Pair", "Epsilon (kJ/mol)", "Sigma (nm)"))
    print("|" + "-"*9 + "+" + "-"*17 + "+" + "-"*17 + "|")
    for item in parsed_data['lj']:
        pair = item['pair']
        eps = item['epsilon'] if item['epsilon'] is not None else float('nan')
        sig = item['sigma'] if item['sigma'] is not None else float('nan')
        pair_str = f"{pair[0]} - {pair[1]}"
        print("| {:<9} | {:>15.7g} | {:>15.7g} |".format(pair_str, eps, sig))
    print("-" * 54)

# ------------------ Analysis Functions ------------------

def error_rdf_analysis(show_plots=True):
    """
    Performs the Error vs. Iteration and RDF descriptor analysis.
    Requires the FOX package.
    """
    import numpy as np
    from FOX import from_hdf5
    from FOX.recipes import overlay_descriptor, plot_descriptor

    hdf5_file = 'armc.hdf5'
    
    # --- Error vs. Iteration Plot ---
    err = from_hdf5(hdf5_file, 'aux_error')
    acceptance = from_hdf5(hdf5_file, 'acceptance')
    accerr = err[acceptance.values].loc[:, "rdf.0"]
    
    plt.figure()
    plt.plot(accerr, label="rdf.0")
    plt.vlines(accerr.idxmin(), accerr.max()+0.5, accerr.min()-0.5, colors='k', linestyles='dashed')
    plt.legend()
    plt.xlabel('ARMC iteration')
    plt.ylabel('Auxiliary Error')
    
    # --- RDF Plot ---
    _dct = overlay_descriptor(hdf5_file, name='rdf')
    keys_list = list(_dct.keys())
    
    print("Available RDF keys (number -> pair):")
    for i, key in enumerate(keys_list, start=1):
        print(f"  {i} -> {key}")
    
    user_input = input(
        "\nEnter the numbers of the RDF keys to plot, separated by spaces.\n"
        "Or press Enter to plot ALL:\n"
    )
    
    if user_input.strip():
        try:
            desired_indices = [int(x) - 1 for x in user_input.split()]
            valid_indices = [idx for idx in desired_indices if 0 <= idx < len(keys_list)]
            selected_keys = [keys_list[idx] for idx in valid_indices]
            dct = {k: _dct[k] for k in selected_keys}
            if not dct:
                print("Warning: No valid selections were found. Plotting nothing.")
        except ValueError:
            print("Error: Invalid input. Plotting nothing.")
            dct = {}
    else:
        dct = _dct
    
    plot_descriptor(
        dct,
        ylim=(-0.5, 14.5),
        yticks=np.arange(0, 14.5, step=2),
        xlabel='Interatomic Distance r (Å)',
        fontsize=20
    )
    
    if show_plots:
        plt.show()

def best_parameters_analysis():
    """
    Finds the move with the minimum error in armc.log, locates the corresponding folder,
    and parses the CP2K input file to print the best parameters.
    """
    folder_str, min_error = get_best_folder("armc.log")
    if folder_str is None:
        print("Could not determine best folder from armc.log.")
        return
    print(f"Minimum error found: {min_error:.4f}")
    print(f"Corresponding folder: md_settings.{folder_str}")
    
    folder_name = f"md_settings.{folder_str}"
    cp2k_in_path = os.path.join(folder_name, f"{folder_name}.in")
    if not os.path.isfile(cp2k_in_path):
        print(f"Error: {cp2k_in_path} not found!")
        return
    parsed_data = parse_cp2k_in(cp2k_in_path)
    print_parameters_table(parsed_data)

def md_analysis(show_plots=True):
    """
    Performs MD analysis (cell lengths, temperature, used time) using the folder corresponding
    to the minimum error move.
    """
    folder_str, _ = get_best_folder("armc.log")
    if folder_str is None:
        print("Could not determine best folder from armc.log.")
        return

    file_md = f"md_settings.{folder_str}/md_settings.{folder_str}.out"
    file_ener = f"md_settings.{folder_str}/cp2k-1.ener"
    
    times_md, cell_x, cell_y, cell_z = [], [], [], []
    
    try:
        with open(file_md, "r") as file:
            current_time = None
            for line in file:
                if "MD| Time [fs]" in line and "Conserved quantity" not in line:
                    parts = line.strip().split()
                    current_time = float(parts[-1])
                elif "MD| Cell lengths [ang]" in line and "Average" not in line:
                    parts = line.strip().split()
                    lx, ly, lz = float(parts[-3]), float(parts[-2]), float(parts[-1])
                    if current_time is not None:
                        times_md.append(current_time)
                        cell_x.append(lx)
                        cell_y.append(ly)
                        cell_z.append(lz)
    except FileNotFoundError:
        print(f"Error: The file '{file_md}' does not exist.")
        return
    
    times_ener, temps, used_time = [], [], []
    
    try:
        with open(file_ener, "r") as file:
            for line in file:
                if line.strip() and not line.startswith("#"):
                    parts = line.strip().split()
                    if len(parts) >= 7:
                        times_ener.append(float(parts[1]))
                        temps.append(float(parts[3]))
                        used_time.append(float(parts[6]))
    except FileNotFoundError:
        print(f"Error: The file '{file_ener}' does not exist.")
        return
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
    
    # Plot cell lengths
    ax1.plot(times_md, cell_x, 'r-', label='Cell Length x (Ang)')
    ax1.plot(times_md, cell_y, 'b-', label='Cell Length y (Ang)')
    ax1.plot(times_md, cell_z, 'g-', label='Cell Length z (Ang)')
    ax1.set_xlabel("Time [fs]")
    ax1.set_ylabel("Cell Length [Ang]")
    ax1.set_title(f"Cell Lengths vs. Time (md_settings.{folder_str})")
    ax1.legend()
    ax1.grid(True)
    
    # Plot Temperature and Used Time with dual y-axes
    ax3 = ax2.twinx()
    ax2.plot(times_ener, temps, 'm-', label='Temperature [K]')
    ax3.plot(times_ener, used_time, 'c-', label='Used Time [s]')
    ax2.set_xlabel("Time [fs]")
    ax2.set_ylabel("Temperature [K]", color='m')
    ax3.set_ylabel("Used Time [s]", color='c')
    ax2.set_title(f"Temperature & Used Time vs. Time (md_settings.{folder_str})")
    ax2.tick_params(axis='y', colors='m')
    ax3.tick_params(axis='y', colors='c')
    ax2.grid(True)
    fig.tight_layout()
    
    if show_plots:
        plt.show()

# ------------------ Main Workflow ------------------

def main():
    print("Choose analysis to perform:")
    print("  1. Error and RDF Analysis")
    print("  2. Print Best Parameters")
    print("  3. MD Analysis")
    print("  4. Run All Analyses")
    choice = input("Enter option number (1, 2, 3, or 4): ").strip()
    
    if choice == "1":
        error_rdf_analysis()
    elif choice == "2":
        best_parameters_analysis()
    elif choice == "3":
        md_analysis()
    elif choice == "4":
        # Run all analyses but suppress individual plt.show() calls;
        # then display all figures together.
        error_rdf_analysis(show_plots=False)
        best_parameters_analysis()
        md_analysis(show_plots=False)
        plt.show()
    else:
        print("Invalid choice. Exiting.")

if __name__ == "__main__":
    main()

