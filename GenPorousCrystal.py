from ase import Atoms
from ase.build import bulk
from ase.io.lammpsdata import write_lammps_data
import numpy as np
from constants import *

def create_unit_cell(lattice_constant=LATTICE_CONSTANT):
    """
    Create a single unit cell of iron with a face-centered cubic (fcc) structure.

    Parameters:
    lattice_constant (float): The lattice constant of the iron crystal structure. Default is LATTICE_CONSTANT.

    Returns:
    Atoms: An ASE Atoms object representing the unit cell.
    """
    return bulk('Fe', 'fcc', a=lattice_constant, cubic=True)

def calculate_atom_numbers(supercell):
    """
    Calculate the number of each type of atom in the supercell.

    Parameters:
    supercell: The supercell for which to calculate atom numbers.

    Returns:
    num_iron, num_chromium, num_nickel, num_carbon (int): The number of iron, chromium, nickel, and carbon atoms.
    """
    # Define the fraction of each type of atom
    fraction_iron = 0.6875
    fraction_chromium = 0.19
    fraction_nickel = 0.12
    fraction_carbon = 0.0025

    # Calculate the total number of substitutional atoms (iron, chromium, nickel)
    fraction_substitutional_atoms = fraction_iron + fraction_chromium + fraction_nickel
    total_num_substitutional_atoms = len(supercell)

    # Calculate the total number of atoms, including interstitial atoms (carbon)
    total_num_atoms = round(total_num_substitutional_atoms / fraction_substitutional_atoms)  # 0.25% interstitials

    # Calculate the number of each type of atom
    num_iron = round(fraction_iron * total_num_atoms)
    num_chromium = round(fraction_chromium * total_num_atoms)
    num_nickel = round(fraction_nickel * total_num_atoms)
    num_carbon = round(fraction_carbon * total_num_atoms)  # approx. 0.25%

    # Check if the total number of substitutional atoms is correct
    total_atoms_check = num_iron + num_chromium + num_nickel

    if total_num_substitutional_atoms != total_atoms_check:
        # Calculate the rounding error
        rounding_error = total_num_substitutional_atoms - total_atoms_check

        # Print the rounding error
        print(f"Rounding error: {rounding_error}")

        # Adjust the number of iron atoms
        print(f"Adjusting number of iron atoms from {num_iron}")
        num_iron = total_num_substitutional_atoms - num_chromium - num_nickel
        print(f"to {num_iron}")
        

    # Return the number of each type of atom
    return num_iron, num_chromium, num_nickel, num_carbon

def create_supercell(num_cells_x=NX, num_cells_y=NY, num_cells_z=NZ, lattice_constant=LATTICE_CONSTANT, include_carbon=INCLUDE_CARBON):
    """
    This function creates a supercell of a crystal structure with specific atom proportions.

    Parameters:
    num_cells_x (int): The number of cells in the x direction. Default is NX.
    num_cells_y (int): The number of cells in the y direction. Default is NY.
    num_cells_z (int): The number of cells in the z direction. Default is NZ.
    lattice_constant (float): The lattice constant of the crystal structure. Default is LATTICE_CONSTANT.
    include_carbon (bool): A boolean to decide whether to include carbon atoms or not. Default is INCLUDE_CARBON.

    Returns:
    supercell: A supercell of the crystal structure with the specified atom proportions.
    """

    def shuffle_and_replace_atoms(supercell, num_iron, num_chromium, num_nickel):
        """Shuffle the atoms and replace the Fe atoms with the other types."""
        
        # Convert supercell to list and shuffle it
        atoms_list = list(supercell)
        np.random.shuffle(atoms_list)

        # Define the ranges for each atom type
        iron_range = range(num_iron)
        chromium_range = range(num_iron, num_iron + num_chromium)
        nickel_range = range(num_iron + num_chromium, num_iron + num_chromium + num_nickel)

        # Replace atoms with the specified types
        for i in iron_range:
            atoms_list[i].symbol = 'Fe'
        for i in chromium_range:
            atoms_list[i].symbol = 'Cr'
        for i in nickel_range:
            atoms_list[i].symbol = 'Ni'

        # Define the cell dimensions
        cell_dimensions = [num_cells_x*lattice_constant, num_cells_y*lattice_constant, num_cells_z*lattice_constant]

        # Create a new supercell with the shuffled atoms and return it
        supercell_shuffled = Atoms(atoms_list, pbc=True, cell=cell_dimensions)
        
        return supercell_shuffled

    def add_carbon_atoms(supercell, num_carbon):
        """Add the carbon atoms at random interstitial sites."""
        
        # Define the positions for the interstitial sites
        interstitial_positions = [
            np.array([0.25, 0.25, 0.25]), 
            np.array([0.75, 0.75, 0.25]), 
            np.array([0.75, 0.25, 0.75]), 
            np.array([0.25, 0.75, 0.75])
        ]
        
        # Loop over the number of carbon atoms to be added
        for _ in range(num_carbon):
            
            # Choose a random interstitial position
            position = np.array(interstitial_positions[np.random.randint(4)])
            
            # Generate random shifts in x, y, and z directions
            shift_x = np.random.randint(num_cells_x)
            shift_y = np.random.randint(num_cells_y)
            shift_z = np.random.randint(num_cells_z)
            
            # Combine the shifts into a single array
            shift = np.array([shift_x, shift_y, shift_z])
            
            # Add the shift to the position
            position += shift
            
            # Scale the position by the lattice constant
            position *= lattice_constant
            
            # Create a new carbon atom at the calculated position
            C = Atoms('C', positions=[position])
            
            # Add the new carbon atom to the supercell
            supercell.extend(C)
        
        # Return the updated supercell
        return supercell

    # Main process of create_supercell

    # Step 1: Create the unit cell with the given lattice constant
    unit_cell = create_unit_cell(lattice_constant)

    # Step 2: Repeat the unit cell in the x, y, and z directions to create the supercell
    supercell = unit_cell.repeat((num_cells_x, num_cells_y, num_cells_z))

    # Step 3: Calculate the number of each type of atom in the supercell
    num_iron, num_chromium, num_nickel, num_carbon = calculate_atom_numbers(supercell)

    # Step 4: Shuffle the atoms in the supercell and replace some of them with other types
    supercell = shuffle_and_replace_atoms(supercell, num_iron, num_chromium, num_nickel)

    # Step 5: Add carbon atoms to the supercell at random interstitial sites
    if include_carbon:
        supercell = add_carbon_atoms(supercell, num_carbon)

    # Return the final supercell
    return supercell

def calculate_delta_distribution_integer_pore_sizes(supercell, porosity=0.00):
    """
    This function calculates the possible integer number of pores and the corresponding number of atoms per pore.
    It also returns the total number of atoms to remove, the discrete porosity, and the porosity error.

    Parameters:
    supercell (list): The supercell to calculate from.
    porosity (float): The porosity of the supercell. Defaults to 0.00.

    Returns:
    tuple: A tuple of five elements. The first two are numpy.ndarrays representing the possible integer number of pores 
        and the corresponding number of atoms per pore. The third is the total number of atoms to remove. 
        The fourth is the discrete porosity. The fifth is the porosity error.
    """

    # Calculate the total number of atoms in the supercell
    total_num_atoms = len(supercell)

    # Calculate the total number of atoms to be removed based on the porosity
    int_num_atoms_to_remove = int(np.round(porosity * total_num_atoms))

    # Calculate the discrete porosity and porosity error
    discrete_porosity = int_num_atoms_to_remove / total_num_atoms
    porosity_error = np.abs(porosity - discrete_porosity) / porosity

    # Generate a list of possible integer numbers of atoms per pore
    possible_atoms_per_pore = np.array([i for i in range(1, int_num_atoms_to_remove + 1)])

    # Calculate the corresponding possible number of pores for each possible number of atoms per pore
    possible_pores = int_num_atoms_to_remove / possible_atoms_per_pore

    # Select only those numbers of pores that are integers
    integer_pores = np.array(possible_pores[possible_pores % 1 == 0],dtype=int)

    # Select the corresponding number of atoms per pore
    integer_atoms_per_pore = np.array(possible_atoms_per_pore[possible_pores % 1 == 0],dtype=int)

    # Return the results
    return integer_pores, integer_atoms_per_pore, int_num_atoms_to_remove, discrete_porosity, porosity_error

def get_atoms_to_remove(supercell, central_atom_index, number_of_atoms_to_remove):
    """
    This function identifies the indices of atoms to remove based on their proximity to a central atom.

    Parameters:
    supercell (ase.Atoms): The supercell to calculate from.
    central_atom_index (int): The index of the central atom.
    cutoff_radius (float): The radius within which to consider atoms for removal.
    number_of_atoms_to_remove (int): The number of atoms to remove.

    Returns:
    list: A list of indices of atoms to remove.
    """

    # Calculate the distances from the central atom to each of the supercell atoms.
    distances = supercell.get_distances(central_atom_index, np.arange(len(supercell)), mic=True)

    # Sort the indices of the supercell atoms based on their distances from the central atom.
    sorted_indices = np.argsort(distances)

    # Use the sorted indices to get a list of the first number of atoms per pore nearby atoms, sorted by distance from the central atom.
    sorted_pore_atom_indices = sorted_indices[:number_of_atoms_to_remove]

    # Return the sorted indices of pore atoms to be removed
    return sorted_pore_atom_indices

def select_pore_configuration(supercell, porosity, choice_index):
    """
    This function selects a specific pore configuration based on a choice index.

    Parameters:
    supercell (ase.Atoms): The supercell to calculate from.
    porosity (float): The desired porosity.
    choice_index (int): The index of the pore configuration to select.

    Returns:
    tuple: A tuple of five elements. The first is the number of pores in the selected configuration. 
        The second is the number of atoms per pore in the selected configuration. 
        The third is the total number of atoms to remove. 
        The fourth is the discrete porosity. 
        The fifth is the porosity error.
    """

    # Calculate the possible pore configurations
    number_of_pores, number_of_atoms_per_pore, total_number_of_atoms_to_remove, discrete_porosity, porosity_error = calculate_delta_distribution_integer_pore_sizes(supercell, porosity)

    # Select the pore configuration based on the choice index
    selected_number_of_pores = number_of_pores[choice_index]
    selected_number_of_atoms_per_pore = number_of_atoms_per_pore[choice_index]

    # Return the selected pore configuration and the other calculated values
    return selected_number_of_pores, selected_number_of_atoms_per_pore, total_number_of_atoms_to_remove, discrete_porosity, porosity_error


def create_porosity_delta_distribution_integer_pore_sizes(supercell, porosity=POROSITY, choice_index=PORE_CONFIGURATION_INDEX):
    """
    This function creates porosity in a supercell by randomly removing atoms.

    Parameters:
    supercell (ase.Atoms): The supercell to modify.
    porosity (float): The desired porosity. Defaults to POROSITY.
    choice_index (int): The index of the pore configuration to select. Defaults to PORE_CONFIGURATION_INDEX.

    Returns:
    tuple: A tuple of five elements. The first is the modified supercell. The second is the ratio of the total number of atoms removed to the total number of atoms.
           The third is the total number of atoms to remove. The fourth is the discrete porosity. The fifth is the porosity error.
    """

    # Select the pore configuration based on the porosity and choice index
    number_of_pores, number_of_atoms_per_pore, total_number_of_atoms_to_remove, discrete_porosity, porosity_error = select_pore_configuration(supercell, porosity, choice_index)

    # Initial total number of atoms
    total_number_of_atoms = len(supercell)
    
    # Initialize the total number of atoms removed
    total_number_of_atoms_removed = 0

    # Create pores
    for _ in range(number_of_pores):

        # Update the number of atoms in supercell
        new_total_number_of_atoms = len(supercell)

        # Randomly select an atom to remove
        atom_to_remove = np.random.randint(new_total_number_of_atoms)

        if number_of_atoms_per_pore > 1:
            # Get the indices of the atoms to remove including the central atom
            indices_to_remove = get_atoms_to_remove(supercell, atom_to_remove, number_of_atoms_per_pore)
        else:
            indices_to_remove = np.array([atom_to_remove])

        
        # Sort the indices in descending order to avoid index shifting during deletion
        indices_to_remove = np.sort(indices_to_remove)[::-1]

        # Remove the atoms
        for atom_i in indices_to_remove:
            supercell.pop(atom_i)

        # Update the total number of atoms removed
        total_number_of_atoms_removed += len(indices_to_remove)

    # Return the modified supercell and the ratio of the total number of atoms removed to the total number of atoms
    return supercell, total_number_of_atoms_removed / total_number_of_atoms, total_number_of_atoms_to_remove, discrete_porosity, porosity_error


def recalculate_num_fraction(supercell):
    """
    This function re-calculates the number fraction of each atom type in the porous supercell.

    Parameters:
    supercell (ase.Atoms): The porous supercell.

    Returns:
    None
    """

    # Get the total number of atoms in the supercell
    total_atoms = len(supercell)

    # Initialize a dictionary to store the number of each atom type
    atom_counts = {'Fe': 0, 'Cr': 0, 'Ni': 0, 'C': 0}

    # Count the number of each atom type
    for atom in supercell:
        atom_counts[atom.symbol] += 1

    # Calculate the number fraction of each atom type and print the results
    for atom_type, count in atom_counts.items():
        num_fraction = count / total_atoms
        print(f'Porous number fraction of {atom_type}: {num_fraction}')

        # Write the results to a text file
        with open('porous_cell_num_fractions.txt', 'a') as file:
            file.write(f'Porous number fraction of {atom_type}: {num_fraction}\n')

def main():
    """
    This function creates a porous supercell and writes it to a LAMMPS data file.
    """

    # Create the supercell
    supercell = create_supercell(NX, NY, NZ)

    if POROSITY == 0.0:
        # Write the supercell to a LAMMPS data file
        if INCLUDE_CARBON:
            write_lammps_data(f'{NX}by{NY}by{NZ}_supercell_0_pores_0P00_porosity.data', supercell,specorder=['Fe','Cr','Ni','C'],atom_style='charge')
        else:
            write_lammps_data(f'{NX}by{NY}by{NZ}_supercell_0_pores_0P00_porosity.data', supercell,specorder=['Fe','Cr','Ni'],atom_style='charge')
        return
    
    # Calculate the possible pore configurations
    possible_pore_numbers, _, _, _, _ = calculate_delta_distribution_integer_pore_sizes(supercell, POROSITY)

    # Print the possible pore number
    print("Possible Number of Pores = ", possible_pore_numbers)

    # Print the number of possible pore configurations
    print("Number of Possibilities for pore sizes/number = ", len(possible_pore_numbers))

    # Print the selected number of pores
    print("Selected Number of Pores = ", possible_pore_numbers[PORE_CONFIGURATION_INDEX])

    # Create porosity in the supercell and get the corrected porosity, error, total number of atoms to remove, and number of atoms per pore
    porous_supercell, created_porosity, total_number_of_atoms_to_remove, discrete_porosity, porosity_error = create_porosity_delta_distribution_integer_pore_sizes(supercell, POROSITY, PORE_CONFIGURATION_INDEX)

    # Print the calculated values
    print("Corrected Discrete Porosity = ", discrete_porosity)
    print("Corrected Discrete Porosity error % = ", porosity_error * 100)
    print("Total number of atoms to remove = ", total_number_of_atoms_to_remove)

    # Print the created porosity
    print("Created Porosity = ", created_porosity)

    recalculate_num_fraction(porous_supercell)

    # Write the porous supercell to a LAMMPS data file
    if INCLUDE_CARBON:
        write_lammps_data(f'{NX}by{NY}by{NZ}_supercell_{possible_pore_numbers[PORE_CONFIGURATION_INDEX]}_pores_0P{int(str(round(created_porosity, 2) % 1)[2:])}_porosity.data', porous_supercell,specorder=['Fe','Cr','Ni','C'],atom_style='charge')
    else:
        write_lammps_data(f'{NX}by{NY}by{NZ}_supercell_{possible_pore_numbers[PORE_CONFIGURATION_INDEX]}_pores_0P{int(str(round(created_porosity, 2) % 1)[2:])}_porosity.data', porous_supercell,specorder=['Fe','Cr','Ni'],atom_style='charge')

    return

if __name__ == "__main__":
    main()