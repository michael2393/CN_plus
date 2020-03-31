"""
|=============================================================================|
|                            I N I T I A L I S E                              |
|=============================================================================|
|                                                                             |
| This module contains routines that set the initial grid and structure       |
| of atoms.                                                                   |
|                                                                             |
| Contains                                                                    |
| --------                                                                    |
|     set_up_unit_cell                                                        |
|     populate_cell_with_atoms                                                |
|     initialise_grid_from_cif    (written by Paul Sharp)                     |
|                                                                             |
|-----------------------------------------------------------------------------|
| Michail Theofilatos 05/11/2019                                              |
|=============================================================================|
"""

import ase
import ase.io

# =============================================================================
# =============================================================================
def set_up_unit_cell(axes_points, cell_spacing, angles, rng):
    """
    Initialise the unit cell. Angles are randomly chosen to be in [minimum,maximum] (given as input: angles=minimum,maximum)

    Parameters
    ----------
    axes_points : integer
        List of the number of points along each dimension.
    cell_spacing : float
        Space between consecutive points on each dimension.
    angles : integer
        List of minimum and maximum angle.
    rng : NR_ran
        Random number generator - algorithm from Numerical Recipes 2007

    Returns
    -------
    initial_struct : ase atoms
        The initial atoms object.

    ---------------------------------------------------------------------------
    Michail Theofilatos 30/10/2019
    """

    #cell_spacing = list(map(float,cell_spacing))
    #axes_points = list(map(int,axes_points))

    print(axes_points)

    initial_struct = ase.Atoms(cell=[float(axes_points[0]) * float(rng.real_range(cell_spacing[0], cell_spacing[1])),
                                     float(axes_points[1]) * float(rng.real_range(cell_spacing[0], cell_spacing[1])),
                                     float(axes_points[2]) * float(rng.real_range(cell_spacing[0], cell_spacing[1])),
                                     rng.int_range(angles[0], angles[1]+1),
                                     rng.int_range(angles[0], angles[1]+1),
                                     rng.int_range(angles[0], angles[1]+1)],
                               pbc=[True, True, True])

    return initial_struct


# =============================================================================
def populate_cell_with_atoms(unit_cell, atoms):
    """
    Places the atoms in random positions inside the unit cell.
    The atoms are treated as hard spheres, thus, overlapping atoms are not allowed.

    Parameters
    ----------
    unit_cell : ase atoms
        The initial structure with the unit cell set.
    atoms : atoms list of structure_class
        Each row of the atoms list consists of the symbol, the number of atoms, the charge and the radius.

    Returns
    -------
    unit_cell : ase atoms
        The structure with the atoms randomly placed in the unit cell.

    ---------------------------------------------------------------------------
    Michail Theofilatos 30/10/2019
    """

    for species_data in atoms:
        positions = []
        for atoms in range(0, species_data["numberOfAtoms"]):
            positions.append([0,0,0])
        unit_cell.extend(ase.Atoms(species_data["symbol"] + str(species_data["numberOfAtoms"]),
                                        cell=unit_cell.get_cell(),
                                        scaled_positions=positions,
                                        charges=[species_data["charge"]] * species_data["numberOfAtoms"],
                                        pbc=[True, True, True]))

    return unit_cell



# =============================================================================
def initialise_from_cif(cif_file, atoms_file_data):
    """
    Initialise an atoms object with a specific structure contained in a cif file.

    ASE cannot read in charge data, so we read it ourselves from either the
    cif file or the ".atoms" file.

    Parameters
    ----------
    cif_file : string
        The ".cif" file containing the initial structure.
    atoms_file_data : str, int, float
        The species, number of atoms, and charge for all atoms in the structure
        from the cif file.

    Returns
    -------
    initial_struct : ase atoms
        The initial structure from the cif file, including charges for each atom.

    ---------------------------------------------------------------------------
    Paul Sharp 22/09/2017
    """

    loop_marker = "loop_"
    symbol_marker = "_atom_site_type_symbol"
    oxidation_number_marker = "_atom_site_type_oxidation_number"

    names = []
    charges = []

    # Read cif file into atoms object
    initial_struct = ase.io.read(cif_file)

    with open(cif_file, mode="r") as cif:
        cif_file_data = cif.readlines()

    # Read charges from either cif file or atoms file
    read_from_cif = False
    for line in cif_file_data:
        if oxidation_number_marker in line:
            read_from_cif = True
            break

    if read_from_cif:

        # Remove blank lines and newline characters
        cif_file_data[:] = [x.strip() for x in cif_file_data if not x.startswith("\n")]

        # Find loops in cif file, and the location of oxidation numbers
        loop_indices = [i for i, x in enumerate(cif_file_data) if x == loop_marker]
        oxidation_index = cif_file_data.index(oxidation_number_marker)

        # Find loop with oxidation numbers
        for index in loop_indices:
            if index > oxidation_index:
                end_index = index
                break
        else:
            end_index = len(cif_file_data)

        try:
            start_index = loop_indices[loop_indices.index(end_index) - 1]
        except ValueError:
            start_index = loop_indices[-1]

        oxidation_loop = cif_file_data[start_index:end_index]

        symbol_pos = oxidation_loop.index(symbol_marker) - 1
        oxidation_pos = oxidation_loop.index(oxidation_number_marker) - 1

        # Underscores are not allowed in atom symbol code -- hence separate header from data
        oxidation_loop_data = [x for x in oxidation_loop if "_" not in x]

        for line in oxidation_loop_data:

            names.append(line.split()[symbol_pos].translate(None, '0123456789'))
            charges.append(float(line.split()[oxidation_pos]))

    # If charge data is not included in the cif file, read it from the atoms file
    else:

        for entry in atoms_file_data:
            names.extend([entry.split()[0]] * int(entry.split()[1]))
            charges.extend([float(entry.split()[2])] * int(entry.split()[1]))

    #print(collections.Counter(names))
    #print(collections.Counter(initial_struct.get_chemical_symbols()))

    # Set charges from list if list matches atoms object
    #if (collections.Counter(names) != collections.Counter(initial_struct.get_chemical_symbols())):
    #    sys.exit("ERROR in initialise.initialise_from_cif() -- the list of atoms with charges does not match up with the atoms in the cif file.")

    initial_struct.set_initial_charges(charges)

    return initial_struct

