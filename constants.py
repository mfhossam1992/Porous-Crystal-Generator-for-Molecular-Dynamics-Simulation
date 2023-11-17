def read_input(filename="input.txt"):
    """
    This function reads an input file and returns a dictionary of constants.

    Parameters:
    filename (str): The name of the input file.

    Returns:
    constants (dict): A dictionary of constants.
    """

    constants = {}

    with open(filename, 'r') as file:
        for line in file:
            # Skip comments and empty lines
            if line[0] == '#' or line.strip() == '':
                continue
            
            name, value = line.strip().split('=')
            # Remove leading and trailing whitespace
            name = name.strip()
            value = value.strip()

            # Convert the value to the appropriate type
            if value.lower() == "true" or value.lower() == "false":
                constants[name] = value.strip().lower() == 'true'
            else:
                constants[name] = float(value) if '.' in value else int(value)
    # Return the dictionary of constants
    return constants

# Read the input file
constants = read_input()

# LATTICE_CONSTANT is a constant that represents the lattice constant of the crystal
LATTICE_CONSTANT = constants['LATTICE_CONSTANT']

# NUM_LATTICE_ATOMS is a constant that represents the approximate desired number of substitutional atoms in the lattice
## Note that this number will get adjusted to the nearest smaller number that fits in a cubic FCC supercell'''
NUM_LATTICE_ATOMS = constants['NUM_LATTICE_ATOMS']

# POROSITY is a constant that represents the desired porosity of the supercell
POROSITY = constants['POROSITY']

# INCLUDE_CARBON is a bool to decide whether to include interstital carbon atoms or not
INCLUDE_CARBON = constants['INCLUDE_CARBON']

# PORE_CONFIGURATION_INDEX is an integer that represents the index of the pore configuration to select
PORE_CONFIGURATION_INDEX = constants['PORE_CONFIGURATION_INDEX']
''' 
The pore configurations are sorted in descending order of number of pores. 
Index 0 is the largest number of pores, each of size 1 atom.
Index 1 is the second largest number of pores, and so on.
Index -1 is the smallest number of pores, which is one big pore.
'''

###############################
# End of user-defined constants
###############################

def calculate_num_cells_per_dimension(total_substitutional_atoms=NUM_LATTICE_ATOMS):
    """
    Calculate the number of cells along each dimension (nx, ny, nz) 
    based on the total number of substitutional atoms.

    Parameters:
    total_substitutional_atoms (int): The desired total number of substitutional atoms in the supercell.
                                      Default is NUM_LATTICE_ATOMS.

    Returns:
    tuple: A tuple containing the number of cells along each dimension (nx, ny, nz).
    """
    # Calculate the number of cells along each dimension
    # The cube root ensures that the supercell is cubic in shape
    num_cells = round((total_substitutional_atoms / 4) ** (1 / 3))

    # The number of cells along each dimension is the same
    nx = ny = nz = num_cells

    return nx, ny, nz

# Calculate the number of cells in each dimension
NX,NY,NZ = calculate_num_cells_per_dimension()
