# Porous-Crystal-Generator-for-Molecular-Dynamics-Simulation

## Overview

This Python code repository provides a tool for generating porous crystal structures suitable for molecular dynamics simulations, with a focus on additive manufactured structures. The tool leverages various algorithms and parameters to create realistic and customizable porous structures, allowing researchers to study the mechanical properties, thermal conductivity, and other characteristics of these materials.

## Features
- **Porous Structure Generation:** The code generates 3D porous crystal structures with user-defined parameters such as porosity, crystal lattice dimensions, and pore size distribution.

- **Additive Manufacturing Considerations:** Special attention is given to the specific requirements of additive manufacturing processes, ensuring that the generated structures are suitable for simulation and, if desired, subsequent physical fabrication.

- **Configurability:** Users can easily configure the crystal structure, porosity levels, and other relevant parameters through a simple input text file.

- **Output Formats:** The generated structures can be exported in common file formats used in molecular dynamics simulations, making them compatible with popular simulation software packages. (Now only LAMMPS format is supported)

## Prerequisites

- Python 3.x
- Required Python packages (listed in `requirements.txt`)

## Installation

1. Clone the repository:
```bash
   git clone https://github.com/mfhossam1992/Porous-Crystal-Generator-for-Molecular-Dynamics-Simulation.git
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

## Usage

1. Configure the parameters in the "input.txt" file to customize the generated porous structure.
2. Make sure that "GenPorousCrystal.py", "constants.py", and "input.txt" are in the directory from which you run the code.
3. Run the main script:
```bash
python GenPorousCrystal.py
```
4. The generated structure will be saved in the specified output file format (now only LAMMPS format is supported) in the same run directory.

## Example

For a quick demo, run the following command:

```bash
python GenPorousCrystal.py
```
This will generate a sample porous crystal structure with parameters read from the input.txt file.

## Contributing

We welcome contributions! Please check the contribution guidelines for details on how to contribute to this project.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

Special thanks to contributors who have participated in the development of this tool.
Hossam Farag thanks Muhammad Abdelghany for support in code development.
Hossam Farag also thanks Mahmoud Mahrous for helpful discussions.

## Contact

For questions or issues, please contact hfarag2@illinois.edu.

