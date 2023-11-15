# Porous-Crystal-Generator-for-Molecular-Dynamics-Simulation

## Overview

This Python code repository provides a tool for generating porous crystal structures suitable for molecular dynamics simulations, with a focus on additive manufactured structures. The tool leverages various algorithms and parameters to create realistic and customizable porous structures, allowing researchers to study the mechanical properties, thermal conductivity, and other characteristics of these materials.

## Features
- **Porous Structure Generation:** The code generates 3D porous crystal structures with user-defined parameters such as porosity, crystal lattice dimensions, and pore size distribution.

- **Additive Manufacturing Considerations:** Special attention is given to the specific requirements of additive manufacturing processes, ensuring that the generated structures are suitable for simulation and, if desired, subsequent physical fabrication.

- **Configurability:** Users can easily configure the crystal structure, porosity levels, and other relevant parameters through a user-friendly interface.

- **Output Formats:** The generated structures can be exported in common file formats used in molecular dynamics simulations, making them compatible with popular simulation software packages.

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

1. Configure the parameters in the config.yml file to customize the generated porous structure.
2. Run the main script:
```bash
python GenPorousCrystal.py
```
3. Follow the on-screen instructions to input any additional parameters and generate the porous crystal structure.
4. The generated structure will be saved in the specified output file format (e.g., XYZ, PDB) in the output/ directory.

## Example

For a quick demo, run the following command:

```bash
python GenPorousCrystal.py --demo
```
This will generate a sample porous crystal structure with default parameters.

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

