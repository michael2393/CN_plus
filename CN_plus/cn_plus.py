"""
|=============================================================================|
|                               C N   P L U S                                 |
|=============================================================================|
|                                                                             |
|                                                                             |
| Contains                                                                    |
| --------                                                                    |
|     CN_PLUS                                                                 |
|     Structure_data                                                          |
|                                                                             |
|-----------------------------------------------------------------------------|
| Michail Theofilatos 30/10/2019                                              |
|=============================================================================|
"""

import os
import structure_class
import initialise
import inputs
import rngs
import functions
import sys
#from Tkinter import *   ## python2


# =============================================================================
params = []
Output = ""

# =============================================================================
# =============================================================================
def CN_PLUS(calc_name):

    # Check that an ".input" file exists, stop execution if not
    input_file = calc_name + ".input"
    if (not os.path.isfile(input_file)):
        sys.exit('ERROR - there is no ".input" file')

    # Set parameters from input file
    params, errors = inputs.parse_input(input_file, calc_name)

    # Convert, and ensure all input file parameters are of the correct type
    params, errors = inputs.convert_values(params, errors)
    params, errors = inputs.handle_dependencies(params, errors)

    if len(errors) > 0:
        inputs.report_input_file_errors(errors, input_file, calc_name + ".error")
        sys.exit("Terminating execution")

    with open(params["output_file"]["value"], mode="a") as output:
        output.write("\n")
        output.write("##############################################################################\n")
        output.write("#                                  CN PLUS                                   #\n")
        output.write("##############################################################################\n")
        output.write("\n")

        output.write("Summary of Inputs\n")
        output.write("\n")
        output.write("-" * 75 + "\n")
        for keyword in params:
            if params[keyword]["specified"]:
                output.write("%s : %s\n" % (keyword, params[keyword]["value"]))
        output.write("-" * 75 + "\n")
        output.write("\n")



        dummy_energy = 0.0
        dummy_volume = 0.0
        dummy_potentials = []
        dummy_derivs = []
        dirpath = os.getcwd()
        print(dirpath)

        # Set GULP command, input keywords and options
        os.environ["GULP_LIB"] = ""
        if params["calculator_cores"]["value"] > 1:
            os.environ["ASE_GULP_COMMAND"] = "timeout " + str(
                params["calculator_time_limit"]["value"]) + " mpirun -np " + str(
                params["calculator_cores"]["value"]) + " " + dirpath + \
                                             params["gulp_executable"]["value"] + " < PREFIX.gin > PREFIX.got"
        else:
            os.environ["ASE_GULP_COMMAND"] = "timeout " + str(
                params["calculator_time_limit"]["value"]) + " " + dirpath + \
                                             params["gulp_executable"]["value"] + " < PREFIX.gin > PREFIX.got"


        rng = rngs.NR_Ran(params["random_seed"]["value"])
        rng.warm_up(params["rng_warm_up"]["value"])


        # =====================================================================
        # Create a 'structure' object and initialize everything...
        s = structure_class.structure_class(params)

        # Check that a ".atoms" file exists, stop if not
        if not os.path.isfile(params["atoms_file"]["value"]):
            sys.exit('ERROR - the ".atoms" file "%s" does not exist.'
                     % (params["atoms_file"]["value"]))

        # Read species, number, charge and radius of ions from ".atoms" file
        with open(params["atoms_file"]["value"], mode="r") as atoms_file:
            for line in atoms_file:
                s.atoms_data.append(line)

        # Initialize atoms data. Store symbol, number, charge and radius in 'atoms' list
        s.initializeAtomsData()

        if params["initial_structure_file"]["specified"]:
            output.write('The initial structure will be read from the file: "%s".\n'
                         % (params["initial_structure_file"]["value"]))

            s.structure = initialise.initialise_from_cif(params["initial_structure_file"]["value"], s.atoms_data)
        else:
            output.write("The atoms will be initialised on a grid.\n")

            # Set up atoms object and grid.
            initial_cell = initialise.set_up_unit_cell(params["grid_points"]["value"], params["grid_spacing"]["value"], params["angles"]["value"], rng)
            s.structure = initialise.populate_cell_with_atoms(initial_cell, s.distinctAtoms)


        current_structure = Structure_data(s.structure, 0, dummy_energy, dummy_volume, dummy_potentials, dummy_derivs)
        current_structure.write_cif("current.cif")



        # Functions

        # optimization - Performs the operations specified in the .input file for a given number of times (phases)
        if params["function"]["value"] == "optimization":
            operations = [params["operation_" + str(i+1)]["value"] for i in range(0, params["operations"]["value"])]
            procedure_parameters = [params["operations"]["value"],
                                    params["phases"]["value"],
                                    operations]
            functions.optimization(s, procedure_parameters, output)

        # compare_basin_hopping_axes - runs the standard basin hopping with and without the axes algorithm as
        # an intermediate step. Plots figures that show the average time to reach specific energy levels.
        if params["function"]["value"] == "compare_basin_hopping_axes":
            procedure_parameters = [-1, params["phases"]["value"], ["Axes", "Basin_hopping"]]
            functions.compare_basin_hopping_axes(s, params["global_minimum_energy"]["value"], procedure_parameters, output)

        # read_results_from_file - Opens the results (in .pkl) obtained in compare_basin_hopping_axes function
        # and plots the results.
        if params["function"]["value"] == "read_results_from_file":
            functions.plot_data(params, params["pkl_file"]["value"], "Axes - Relaxation")

    return None




# =============================================================================
# =============================================================================
class Structure_data(object):
    """
    Stores data referring to particular structures.

    ---------------------------------------------------------------------------
    Paul Sharp 04/08/2017
    """

    # =========================================================================
    def __init__(self, atoms, index, energy, volume, potentials, derivs):
        """
        Initialise the Structure data.

        Parameters
        ----------
        atoms : ase atoms
            The atoms object containing the structure.
        index : int
            The number of basin hopping moves to get to this structure.
        energy : float
            The energy of the structure.
        volume : float
            The volume of the structure.
        potentials : float
            Site potential for each atom in the current structure.
        derivs : float
            Resolved derivatives of the site potentials for each atom in the current structure.

        Returns
        -------
        None

        -----------------------------------------------------------------------
        Paul Sharp 04/08/2017
        """

        self.atoms = atoms
        self.index = index
        self.energy = energy
        self.volume = volume
        self.potentials = potentials
        self.derivs = derivs

    # =========================================================================
    def write_cif(self, filename):
        """
        Write the structure to a cif file.

        Parameters
        ----------
        filename : str
            The name of the cif file.

        Returns
        -------
        None

        -----------------------------------------------------------------------
        Paul Sharp 03/08/2017
        """

        self.atoms.write(filename, format="cif")

