"""
|=====================================================================================|
|                                       CN PLUS                                       |
|=====================================================================================|
|                                                                                     |
| Welcome to CN PLUS.                                                                 |
| In order to use the code, type:                                                     |
|                                                                                     |
|         python cn_plus file.input                                                   |
|                                                                                     |
|-------------------------------------------------------------------------------------|
| Michail Theofilatos 30/10/2019                                                      |
| I would like to thank Paul Sharp for his help. Many parts of this project	          |
| were developed by him                                                               |
| MODIFIED - Parts of main, cn_plus, inputs, gulp_calc                                |
| UNMODIFIED - Methods/classes: Structure_data, symmetrise_atoms, initialise_from_cif |
|              Files:   rrngs.py                                                      |
|								                                                      |
|=====================================================================================|
"""

import argparse
import ase.io
import os
import sys
import inputs
import cn_plus
try:
    import spglib
except ImportError:
    from pyspglib import spglib



# =============================================================================
# =============================================================================
def main(args=None):

    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser(description='CN PLUS',
                                     epilog='This code is developed by Michail Theofilatos.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('seed', nargs='?', default=None,
                        help='The seed name for this calculation, which forms the name of the ".input" and ".atoms" files.')
    parser.add_argument('-i', '--input', action="store_true", default=False,
                        help='Print all options for the ".input" file with a description of each option.')
    parser.add_argument('-p', '--parse', nargs='+', metavar="<input file>",
                        help='Parse the given input file, report any errors and exit.')
    parser.add_argument('-w', '--write', nargs='?', const="defaults.input",
                        metavar="<input file>",
                        help='Write an input file that includes all keywords with their default values to the given file and exit.')
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s 1.2')

    opts = parser.parse_args(args)

    if opts.input:

        params = inputs.initialise_default_params("")
        print("")

        for key in sorted(params.keys()):
            print("%-30s %s" % (key, params[key]["description"]))
        print("")

    if opts.parse is not None:

        for input_file in opts.parse:

            # Check if file exists
            if not os.path.isfile(input_file):
                print('WARNING - the file "%s" does not exist' % (input_file))
                continue

            calc_name = input_file.split(".")[0]
            params, errors = inputs.parse_input(input_file, calc_name)

            # Convert, and ensure all input file parameters are of the correct type
            params, errors = inputs.convert_values(params, errors)
            params, errors = inputs.handle_dependencies(params, errors)

            if len(errors) > 0:
                print('%i errors were found in the input file: "%s"\n' % (len(errors), input_file))
                print('Please refer to the error log "%s".\n' % (calc_name + ".error"))
                inputs.report_input_file_errors(errors, input_file, calc_name + ".error")
            else:
                print('The file "%s" is a valid input file.' % (input_file))
                print("")
                print("\tPlease note that keywords, options and settings for GULP/VASP will be verified in GULP/VASP calculations.")
                print("\tAlso, when the initial structure is composed, the validity of any supplied swap groups will be verified at that point.")
                print("")


    if opts.write is not None:

        # Check file exists, call write routine if not
        if not os.path.isfile(opts.write):
            inputs.write_defaults_to_file(opts.write)
            print("")
            print('Default values written to input file "%s"' % (opts.write))
            print("")

        else:
            print("")
            print('ERROR - the file "%s" already exists' % (opts.write))
            print("")



    if opts.seed is not None:
        calc_seed = opts.seed.split(".")[0]
        cn_plus.CN_PLUS(calc_seed)

    return None



def symmetrise_atoms(atoms):
    """
    Use spglib in order to find a standardised unit cell from an atoms object.

    Parameters
    ----------
    atoms : ase atoms
        The atoms object containing the structure for which we wish to find symmetries.

    Returns
    -------
    symmetrised_atoms : ase atoms
        The atoms object after finding symmetries.

    ---------------------------------------------------------------------------
    Paul Sharp 27/10/2017
    """

    # Call spglib on the atoms object
    lattice, scaled_pos, atomic_numbers = spglib.standardize_cell(atoms)

    # Record results in a new atoms object
    try:
        symmetrised_atoms = ase.Atoms(cell=lattice,
                                      scaled_positions=scaled_pos,
                                      numbers=atomic_numbers,
                                      charges=atoms.get_initial_charges(),
                                      pbc=[True, True, True])
    except:
        symmetrised_atoms = ase.Atoms(cell=lattice,
                                      scaled_positions=scaled_pos,
                                      numbers=atomic_numbers,
                                      pbc=[True, True, True])

    return symmetrised_atoms


# =============================================================================
# =============================================================================


if __name__ == "__main__":
    main()
