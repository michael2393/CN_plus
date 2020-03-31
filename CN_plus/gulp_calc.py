"""
|=============================================================================|
|                             G U L P   C A L C                               |
|=============================================================================|
|                                                                             |
| This module contains routines that run custom GULP calculations, and check  |
| the status of GULP calculations.                                            |
|                                                                             |
| Contains                                                                    |
| --------                                                                    |
|     gulp_relaxation                                                         |
|     run_gulp                                                                |
|     read_from_restart_file                                                  |
|     create_input_file_from_restart_file                                     |
|     execute_gulp_command_or_script                                          |
|     read_energy                                                             |
|     read_gnorm                                                              |
|     check_convergence                                                       |
|     check_float                                                             |
|     check_timed_out                                                         |
|     read_outcome                                                            |
|     update_outcomes                                                         |
|     read_potentials                                                         |
|     strip_vacancies                                                         |
|                                                                             |
|-----------------------------------------------------------------------------|
| Paul Sharp 08/03/2018                                                       |
|=============================================================================|
"""

try:
    from ase.calculators.gulp import GULP
except ImportError:
    from ase.calculators.gulp import Gulp as GULP

import math
import numpy as np
import os
import copy
import sys
import re
import subprocess
import time
import shutil

# =============================================================================
# =============================================================================
def gulp_relaxation(s, operation, gulp_shells=[]):
    gulp_start_time = time.time()
    structure = s.structure.copy()
    result = energy = gnorm = ""
    failed = ["failed", "timed out", "unconverged"]
    timeout_outcome = "CPU limit has been exceeded - restart optimisation "
    files = {
        "name": "temp",
        "input" : "temp.gin",
        "output" : "temp.got",
        "restart" : "temp.res"
    }


    # Set calculator
    if operation == "relax_structure":
        s.relax_structure_options.append("dump " + files["restart"])
        calc = (GULP(label=files["name"], keywords=s.relax_structure_keywords,
                 options=s.relax_structure_options, shel=gulp_shells,
                 library=s.gulp_library))
    elif operation == "relax_unit_cell":
        s.relax_unit_cell_options.append("dump " + files["restart"])
        calc = (GULP(label=files["name"], keywords=s.relax_unit_cell_keywords,
                     options=s.relax_unit_cell_options, shel=gulp_shells,
                     library=s.gulp_library))
    structure.set_calculator(calc)

    # Run calculation
    try:
        energy = structure.get_potential_energy()
    except:
        pass
    else:
        gnorm = structure.calc.get_Gnorm()


    # Check if failed, converged, timed out
    if not os.path.isfile(files["output"]):
        result = "gulp failure"
    if read_outcome(files["output"]) == timeout_outcome or check_timed_out(files["output"]):
        result = "timed out"
    if (not isinstance(energy, float)) or (not isinstance(gnorm, float)) or (gnorm > s.max_gnorm):
        result = "unconverged"

    # This seems to be the only way to obtain the relaxed structure...
    try:
        structure, energy = read_from_restart_file(structure, energy, files["restart"])
    except ValueError:
        result = "unconverged"

    if result not in failed:
        if operation == "relax_structure" and check_convergence(files["output"], s.relax_structure_keywords):
            result = "converged"
        elif operation == "relax_unit_cell" and check_convergence(files["output"], s.relax_unit_cell_keywords):
            result = "converged"
        else:
            result = "unconverged"
    if result == "timed out" or result == "gulp failure":
        if os.path.isfile(files["output"]):
            result += " - " + read_outcome(files["output"])

    calc_time = time.time() - gulp_start_time

    return structure, energy, result, calc_time





# =============================================================================
def read_from_restart_file(structure, energy, gulp_res_file):
    """
    Read unit cell, atomic positions and energy from a GULP ".res" file, where
    they are quoted to greater precision than the output file (and hence the
    ASE atoms object if available).

    For the GULP calculator in ASE 3.14-, the unit cell and atomic positions
    (if in fractional coordinates) are not recorded after a calculation and so
    must be obtained from a restart file (or the output).

    Parameters
    ----------
    structure : ASE atoms
        The structure used in the GULP calculation.
    energy : float
        Energy of the structure in eV from the calculation.
    gulp_res_file : string
        The restart file written from this GULP calculation.

    Returns
    -------
    structure : ASE atoms
        The structure used in the GULP calculation with unit cell and atomic
        positions taken from the restart file.
    energy : float
        A high-precision value of the final energy of the calculated structure,
        converted to units of eV/atom.

    ---------------------------------------------------------------------------
    Paul Sharp 15/01/2018
    """

    if os.path.isfile(gulp_res_file):

        with open(gulp_res_file, mode="r") as f:
            res_file = f.readlines()

        for i, line in enumerate(res_file):

            # This matches only if the line consists of "cell" and nothing else
            # Hence, we avoid a match on the "cellonly" keyword
            if re.search("^cell$", line):
                cell_line = res_file[i + 1].split()
                structure.set_cell(np.array([float(cell_line[0]), float(cell_line[1]),
                                             float(cell_line[2]), float(cell_line[3]),
                                             float(cell_line[4]), float(cell_line[5])]))

            # This matches only if the line consists of "vectors" and nothing
            # else. Hence, we avoid a match on the "eigenvectors" keyword
            if re.search("^vectors$", line):

                cell_vectors = []
                for j in range(1, 4):
                    vector_line = res_file[i + j].split()
                    cell_vectors.append([float(vector_line[0]),
                                         float(vector_line[1]),
                                         float(vector_line[2])])

                cell_vectors = np.array(cell_vectors)
                structure.set_cell(cell_vectors)

            # This may match on either the coordinate type as desired, or the
            # "cartesian" keyword. If it matches on the keyword it will break
            # the loop immediately as neither "core" nor "shell" should be on
            # the next line.
            if "cartesian" in line or "fractional" in line:

                coordinate_type = line.split()[0]
                pos = []

                for j in range(i + 1, len(res_file)):

                    if "core" in res_file[j]:

                        atom_pos = []
                        atom_line = res_file[j].split()

                        for line_index in range(2, 5):

                            coord = atom_line[line_index]
                            if "/" in coord:
                                atom_pos.append(float(coord.split("/")[0]) / float(coord.split("/")[1]))
                            else:
                                atom_pos.append(float(coord))

                        pos.append(atom_pos)

                    elif "shell" in res_file[j]:
                        continue

                    else:
                        break

                pos = np.array(pos)

            if "totalenergy" in line:
                try:
                    energy = float(line.split()[1]) / float(len(structure.get_positions()))
                except ValueError:
                    energy = line.split()[1]

                    # Replace positions of atoms in structure, depending on coordinate type
        try:
            coordinate_type
        except NameError:
            print('WARNING in "gulp_calc.read_from_restart_file()" -- no valid coordinate type found in ".res" file')
        else:
            if coordinate_type == "fractional":
                structure.set_scaled_positions(pos)
            elif coordinate_type == "cartesian":
                structure.set_positions(pos)
            else:
                print(
                    'WARNING in "gulp_calc.read_from_restart_file()" -- no valid coordinate type found in ".res" file')

    return structure, energy


# =============================================================================
def check_convergence(gulp_out_file, gulp_keywords):
    """
    Examines an output file of a GULP calculation to determine whether or not a
    structure optimisation was successful.

    Parameters
    ----------
    gulp_out_file : string
        The output file of a GULP calculation
    gulp_keywords : string
        The keywords used for this GULP calculation

    Returns
    -------
    converged : boolean
        Determines whether the GULP calculation has successfully optimised the
        structure

    ---------------------------------------------------------------------------
    Paul Sharp 22/09/2016
    """

    optimisation_keywords = ["opti", "optimise", "grad", "gradient", "fit"]

    if any(keyword in gulp_keywords for keyword in optimisation_keywords):
        optimisation_marker = "**** Optimisation achieved ****"
    else:
        optimisation_marker = "Components of energy"

    converged = False

    with open(gulp_out_file, mode="r") as out_file:

        for line in out_file:
            if optimisation_marker in line:
                converged = True
                break

    return converged


# =============================================================================
def check_float(value):
    representable = True

    try:
        float(value)
    except ValueError:
        representable = False

    return representable


# =============================================================================
def check_timed_out(gulp_out_file):
    finished_marker = "Job Finished"
    terminated_marker = "Program terminated"

    final_lines = subprocess.check_output(["tail", "-2", gulp_out_file])
    timed_out = True
    if (finished_marker in final_lines) or (terminated_marker in final_lines):
        timed_out = False

    return timed_out


# =============================================================================
def read_outcome(gulp_out_file):
    outcome_marker = " **** "
    outcome = ""

    with open(gulp_out_file, mode="r") as out_file:

        for line in out_file:
            if outcome_marker in line:
                outcome = outcome + line.split(outcome_marker)[1][:-5]

    return outcome.strip()
