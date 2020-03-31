# CN_plus
Local search algorithms for the Crystal Structure Prediction problem.

This project is part of the paper published in SEA 2020 (arXiv version: https://arxiv.org/abs/2003.12442).


Arguments
===============================================

python CN_plus -i
Print all options for the ".input" file with a description of each option.

python CN_plus -p path_to_input_file.input
Parse the given input file, report any errors and exit.

python CN_plus -w
Write an input file that includes all keywords with their default values to the given file and exit.


Execution:
===============================================

You can run the code in Python 2.7. It uses the GULP (http://gulp.curtin.edu.au/gulp/) library to perform energy calculations and structure optimization. You can download it here https://gulp.curtin.edu.au/gulp/request.cfm?rel=download.
This project was implemented using GULP 5.1.

After you download GULP, extract it and place it in the same directory as CN_plus. In addition, in the .inputs file, set the keyword "gulp_executable" to point the gulp executable (e.g., gulp_executable=/gulp-5.1/Src/gulp).

python CN_plus Input/SrTiO3/STO.input


Input files:
===============================================

In the folder CN_plus/Input there are two examples of input files and the corresponding potentials (.lib files).

.atoms file: This file defines the type and number of atoms in the unit cell, their charge and their atomic radius.

.lib file: This file specifies the buckingham potential parameters.

.input file: In this file, you can define the initial size and shape of the unit cell, the discretization (for the neighbourhood search), and the GULP input parameters (such as the .lib file, and the gulp options and keywords). You can find information about these in https://gulp.curtin.edu.au/gulp/help/new_help_40_txt.html


Functions:
===============================================

In the .input file, you can choose which function of the code to run, such as optimization, compare_basin_hopping_axes, and open_csv_results.

optimization: This function runs the operations that you specify and their order. The operations that you can run are:
(1) random_structure: Creates a random (valid) structure (number and types of atoms as in .atoms file, and unit cell as defined in the .input file).
(2) axes, swap_all, swap_atoms: Runs the algorithm (greedy) as described in the paper "Crystal Structure Prediction via Oblivious Local Search" (https://arxiv.org/abs/2003.12442), using the corresponding neighbourhood.
(3) relax_structure: Runs GULP optimization in order find the local minimum.
(4) relax_unit_cell: Runs GULP optimization but only changes the unit cell parameters (not the positions of the atoms).
(5) input_structure: Reads a .cif file from the folder indicated by the keyword "initial_structure_folder" in the .inputs file.


Operations:
===============================================

The "operations" keyword specifies the number of operations you want to execute, and the "phases" keyword the number of times to execute them.
The "operation_{i}" keywords defines the operation to run in step i (one of the above).

The progress of the simulation is stored in the file specified by the keyword "output_file" in the .inputs file.



Example:
===============================================

operations = 3
operation_1 = random_structure
operation_2 = axes
operation_3 = relax_structure
phases = 1000
 

This will first generate a random structure, it will run the "axes neighbourhood" algorithm, and then it will use GULP to find the local minimum.
This process will be executed 1000 times.

When it finishes, the results are saved in a csv file in the same directory with the "output_file" (specified in the .inputs file, e.g., /Results/{Date_Time}), and all the structures of each operation in the corresponding folder (e.g.,  /Results/{Date_Time}/random_structure).
The "input_structure" operation reads the file "name_{index}.cif", where index is the current phase of the algorithm (i.e., if you have 10 phases, you need 10 input files).
