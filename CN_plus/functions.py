"""
|=============================================================================|
|                                  FUNCTIONS                                  |
|=============================================================================|
|                                                                             |
|                                                                             |
| Contains                                                                    |
| --------                                                                    |
|     Progress                                                                |
|     optimization                                                            |
|     writeProgressToCSV                                                      |
|     runAxesHeuristic                                                        |
|     relaxStructure                                                          |
|     randomValidStructure                                                    |
|     make_folders                                                            |
|     save_structure                                                          |
|     say                                                                     |
|                                                                             |
|-----------------------------------------------------------------------------|
| Michail Theofilatos 04/11/2019                                              |
|=============================================================================|
"""


try:
    from ase.calculators.gulp import GULP
except ImportError:
    from ase.calculators.gulp import Gulp as GULP

import rngs
from datetime import datetime
import os
import initialise
import gulp_calc
import copy
#import matplotlib.pyplot as plt
#import matplotlib.patches as patches


#from Tkinter import *   ## python2

import heuristics
import cn_plus
import time
import glob
import csv
import statistics
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import sys

# =============================================================================
# =============================================================================

output = ""


class Progress(object):
    """
    Class that stores the progress of the 'optimization' algorithm.

    ---------------------------------------------------------------
    Michail Theofilatos 04/11/2019
    """
    def __init__(self, phase, operation, energy, time, structure, message, numberOfLoops = 0, numberOfSinglePointCalculations = 0):
        self.phase = phase
        self.operation = operation
        self.Energy = energy
        self.Time = time
        self.structure = structure
        self.message = message
        self.numberOfLoops = numberOfLoops
        self.numberOfSinglePointCalculations = numberOfSinglePointCalculations



def optimization(s, procedure_parameters, out):
    """
    This function runs the operation specified in the .input file for 'phases' number of times.

    Parameters
    ----------
    s : structure_class
        Object of structure_class which contains the atoms, the configuration, etc.
    procedure_parameters : [operations, phases, [list of operations]]
        Describes the operations to be executed (e.g., axes, gulp relaxation etc.).
    out : file
        File to write the output text

    Returns
    -------
    s : structure_class

    ---------------------------------------------------------------------------
    Michail Theofilatos 04/11/2019
    """

    global output
    output = out
    progress = []
    initS = copy.copy(s)

    base_path = os.path.join(s.params["output_file"]["value"].rpartition('/')[0], str(datetime.now().strftime("%d_%m_%Y_%H_%M_%S")))
    path = make_folders(base_path, procedure_parameters)


    for phase in range(0, procedure_parameters[1]):
        say("\n\n\n\n Phase " + str(phase) + " \n===========\n")
        s = copy.copy(initS)

        for operation in procedure_parameters[2]:
            initTime = time.time()
            say("\n==== " + operation + " ====")

            # Computes a valid random configuration of atoms
            if operation == "random_structure":
                s = randomValidStructure(s)
                progress.append(Progress(phase = phase,
                                         operation = operation,
                                         energy = s.getEnergy(True)[1],
                                         time = time.time() - initTime,
                                         structure = s.structure,
                                         message = "Found a valid random configuration of atoms"))


            # Performs a gulp calculation in order to find a local minimum
            if operation == "relax_structure" or operation == "relax_unit_cell":
                s, success, result = relaxStructure(s, operation)
                print("Result: " + result + " - energy: " + str(s.getEnergy(True)[1]))
                progress.append(Progress(phase=phase,
                                         operation=operation,
                                         energy=s.getEnergy(False)[1],
                                         time=time.time()-initTime,
                                         structure=s.structure,
                                         message=result))

            # Performs local search using the 'axes' neighbourhood (greedy strategy)
            if operation == "axes":
                s, success, numberOfLoops, numberOfSinglePointCalculations, EnergiesTimes = heuristics.axisHeuristic(copy.copy(s))
                progress.append(Progress(phase=phase,
                                         operation=operation,
                                         energy=s.getEnergy(True)[1],
                                         time=time.time()-initTime,
                                         structure=s.structure,
                                         message="Improved structure: " + str(success),
                                         numberOfLoops=numberOfLoops,
                                         numberOfSinglePointCalculations=numberOfSinglePointCalculations))

            # Performs local search using the 'Swap all' neighbourhood (greedy strategy)
            if operation == "swap_all":
                s, success, numberOfLoops, numberOfSinglePointCalculations = heuristics.swapAllHeuristic(copy.copy(s))
                progress.append(Progress(phase=phase,
                                         operation=operation,
                                         energy=s.getEnergy(True)[1],
                                         time=time.time()-initTime,
                                         structure=s.structure,
                                         message="Improved structure: " + str(success),
                                         numberOfLoops=numberOfLoops,
                                         numberOfSinglePointCalculations=numberOfSinglePointCalculations))

            # Performs local search using the 'Swap atoms' neighbourhood (greedy strategy)
            if operation == "swap_atoms":
                s, success, numberOfLoops, numberOfSinglePointCalculations = heuristics.swapAtomsHeuristic(copy.copy(s))
                progress.append(Progress(phase=phase,
                                         operation=operation,
                                         energy=s.getEnergy(True)[1],
                                         time=time.time()-initTime,
                                         structure=s.structure,
                                         message="Improved structure: " + str(success),
                                         numberOfLoops=numberOfLoops,
                                         numberOfSinglePointCalculations=numberOfSinglePointCalculations))

            # Reads .cif file from Input/Input_Structures
            # In particular, all input files are of the form name_{index}.cif
            # It reads the file with index=phase
            if operation == "input_structure":
                files = [f for f in glob.glob(s.params["initial_structure_folder"]["value"] + "*.cif")]
                structure = filter(lambda k: ("_" + str(phase) + ".") in k, files)
                if(len(structure) > 0):
                    s.structure = initialise.initialise_from_cif(structure[0], s.atoms_data)
                    say("Energy of input structure: " + str(s.getEnergy(True)))
                    say("Read input structure: " + structure[0])
                    progress.append(Progress(phase=phase,
                                             operation=operation,
                                             energy=s.getEnergy(True)[1],
                                             time=time.time() - initTime,
                                             structure=s.structure,
                                             message="Input structure: " + str(structure[0])))
                else:
                    sys.exit("Error: phases are more than the input structures in " + s.params["initial_structure_folder"]["value"])



            save_structure(s, path=path + '/' + operation, name=operation + '_', index=phase)

    # Write results in csv file
    writeProgressToCSV(progress, procedure_parameters, path, name="")




def compare_basin_hopping_axes(initial_s, global_minimum_energy, procedure_parameters, out):
    """
    This function compares the standard basin hopping algorithm with and without the axes neighbourhood
    algorithm as an intermediate step, and saves the results in {output folder}/Compare_Basin_Axes_{date_time}.
    Each iteration consists of the following steps. Computes a valid random configuration of atoms,
    in the first case performs a relaxation in order to find the local minimum, and in the second case
    performs a local search using the 'axes' neighbourhood and then a relaxation. The number of iterations
    is specified by the 'phases' parameter in the .input file.

    Parameters
    ----------
    s : structure_class
        Object of structure_class which contains the atoms, the configuration, etc.
    global_minimum_energy : float
        Each iteration terminates once a configuration of atoms with energy equal or less than
        this parameter is found.
    procedure_parameters : [operations, phases, [list of operations]]
        Describes the operations to be executed (e.g., axes, gulp relaxation etc.).
    out : file
        File to write the output text

    ---------------------------------------------------------------------------
    Michail Theofilatos 10/11/2019
    """

    global output
    output = out

    progress_basin = []
    progress_axes = []
    TimesAndEnergies_axes = []
    TimesAndEnergies_basin = []

    total_loops = 0
    pp_tmp = [2, 0, ["Axes", "Basin_hopping"]]
    base_path = os.path.join(initial_s.params["output_file"]["value"].rpartition('/')[0], str("Compare_Basin_Axes_" + datetime.now().strftime("%d_%m_%Y_%H_%M_%S")))
    path = make_folders(base_path, pp_tmp)

    pp_tmp = [3, 0, ["random_structure", "axes", "relax_structure"]]
    path_axes = make_folders(path + '/Axes', pp_tmp)
    pp_tmp = [2, 0, ["random_structure", "relax_structure"]]
    path_basin = make_folders(path + '/Basin_hopping', pp_tmp)


    # Axes
    for phase in range(0, procedure_parameters[1]):
        print("=== Experiment " + str(phase) + " ===")
        found = False
        loop = 0
        TimesAndEnergies_axes.append([])

        while(not found):
            initTime = time.time()
            print("=== Loop " + str(loop) + " ===")
            s = randomValidStructure(copy.copy(initial_s))

            progress_axes.append(Progress(phase=total_loops,
                                                 operation="random_structure",
                                                 energy=s.getEnergy(True)[1],
                                                 time=time.time() - initTime,
                                                 structure=s.structure,
                                                 message=""))
            TimesAndEnergies_axes[phase].append([s.getEnergy(True)[1], time.time()])
            save_structure(s, path=path_axes + '/random_structure', name="random_structure" + '_', index=total_loops)

            initTime = time.time()
            s, success, numberOfLoops, numberOfSinglePointCalculations, Energies_Times = heuristics.axisHeuristic(copy.copy(s))
            progress_axes.append(Progress(phase=total_loops,
                                                 operation="axes",
                                                 energy=s.getEnergy(True)[1],
                                                 time=time.time() - initTime,
                                                 structure=s.structure,
                                                 numberOfLoops=numberOfLoops,
                                                 message="",
                                                 numberOfSinglePointCalculations=numberOfSinglePointCalculations))
            for p in Energies_Times:
                TimesAndEnergies_axes[phase].append([p[0], p[1]])
            save_structure(s, path=path_axes + '/axes', name="axes" + '_', index=total_loops)


            initTime = time.time()
            s, success, result = relaxStructure(s, "relax_structure")
            print("Result: " + result + " - energy: " + str(s.getEnergy(True)[1]))
            progress_axes.append(Progress(phase=total_loops,
                                     operation="relax_structure",
                                     energy=s.getEnergy(False)[1],
                                     time=time.time() - initTime,
                                     structure=s.structure,
                                     message=result))
            TimesAndEnergies_axes[phase].append([s.getEnergy(True)[1], time.time()])
            save_structure(s, path=path_axes + '/relax_structure', name="relax_structure" + '_', index=total_loops)

            if(s.getEnergy(True)[1]<global_minimum_energy):
                found = True
                print("== FOUND ==")

            loop += 1
            total_loops += 1

    pp_tmp = [3, total_loops, ["random_structure", "axes", "relax_structure"]]
    writeProgressToCSV(progress_axes, pp_tmp, path, name="_axes")
    with open(path + '/Energies_Times_Axes.pkl', 'wb') as f:
        pickle.dump([TimesAndEnergies_axes], f)

    total_loops = 0

    # Basin hopping
    for phase in range(0, procedure_parameters[1]):
        print("=== Experiment " + str(phase) + " ===")
        found = False
        TimesAndEnergies_basin.append([])

        loop = 0

        while (not found):
            print("=== Loop " + str(loop) + " ===")

            initTime = time.time()
            s = randomValidStructure(copy.copy(initial_s))
            progress_basin.append(Progress(phase=total_loops,
                                                 operation="random_structure",
                                                 energy=s.getEnergy(False)[1],
                                                 time=time.time() - initTime,
                                                 structure=s.structure,
                                                  message=""))
            TimesAndEnergies_basin[phase].append([s.getEnergy(True)[1], time.time()])
            save_structure(s, path=path_basin + '/random_structure', name="random_structure" + '_', index=total_loops)

            initTime = time.time()
            s, success, result = relaxStructure(s, "relax_structure")
            print("Result: " + result + " - energy: " + str(s.getEnergy(True)[1]))
            progress_basin.append(Progress(phase=total_loops,
                                                 operation="relax_structure",
                                                 energy=s.getEnergy(False)[1],
                                                 time=time.time() - initTime,
                                                 structure=s.structure,
                                                 message=result))
            TimesAndEnergies_basin[phase].append([s.getEnergy(True)[1], time.time()])
            save_structure(s, path=path_basin + '/relax_structure', name="relax_structure" + '_', index=total_loops)

            if (s.getEnergy(True)[1] < global_minimum_energy):
                found = True
                print("== FOUND ==")
            loop += 1
            total_loops += 1


    pp_tmp = [2, total_loops, ["random_structure", "relax_structure"]]
    writeProgressToCSV(progress_basin, pp_tmp, path, name="_basin_hopping")
    with open(path + '/Energies_Times_Basin_Hopaping.pkl', 'wb') as f:
        pickle.dump([TimesAndEnergies_basin], f)

    say("==== AXES ====")
    compute_results(progress_axes, global_minimum_energy)
    say("==== BASIN HOPPING ====")
    compute_results(progress_basin, global_minimum_energy)


    plot_data(s.params, path + '/Energies_Times_Axes.pkl', name="Axes - Relaxation")
    plot_data(s.params, path + '/Energies_Times_Basin_Hopaping.pkl', name="Basin Hopping")




def compute_results(progress_array, global_minimum_energy):
    """
    Computes the results of the comparisson between basin hopping with and without the intermediate local search (axes neighbourhood)

    progress_array: list of Progress objects
        Contains all the details of the execution.
    global_minimum_energy: float
        The energy of the configuration of atoms that corresponds to the global optimum.
    """

    relaxations_tmp = 0
    relaxations = []
    time_axes = []
    time_relaxations = []
    total_time = []
    total_time_tmp = 0
    time_axes_tmp = 0
    time_relaxations_tmp = 0

    for p in progress_array:
        total_time_tmp += p.Time
        if p.operation == "axes":
            time_axes_tmp += p.Time
        if p.operation == "relax_structure":
            time_relaxations_tmp += p.Time
            relaxations_tmp += 1

        if p.Energy <= global_minimum_energy:
            total_time.append(total_time_tmp)
            time_axes.append(time_axes_tmp)
            time_relaxations.append(time_relaxations_tmp)
            relaxations.append(relaxations_tmp)
            total_time_tmp = 0
            time_relaxations_tmp = 0
            time_axes_tmp = 0
            relaxations_tmp = 0

    print_results(relaxations, "(relaxations)")
    print_results(time_axes, "(axes time)")
    print_results(time_relaxations, "(relaxations time)")
    print_results(total_time, "(total time)")




def print_results(array, name):
    say("\nMean " + name + ": " + str(statistics.mean([i for i in array])))
    say("\nMedian " + name + ": " + str(statistics.median([i for i in array])))
    say("\nstdev " + name + ": " + str(statistics.stdev([i for i in array])))
    say("\nVariance " + name + ": " + str(statistics.variance([i for i in array])))
    say("\n")




def plot_data(params, path, name):
    """
    Creates two 'whisker' plots that show the time needed for the energy to reach specific levels.

    path: string
        Path to the .pkl files produced by 'compare_basin_hopping_axes' function.
    name: string
        Title of the plot.
    """

    with open(path, 'rb') as f1:
        TimesAndEnergies = pickle.load(f1)[0]
        ranges1 = params["whiskers_fine"]["value"]      #[-30.0, -30.4, -30.8, -31.2, -31.4, -31.6, -31.69]
        ranges2 = params["whiskers_coarse"]["value"]    #[-15, -20, -25, -30]

        data1 = [[] for x in ranges1]


        for experiment in TimesAndEnergies:  # For each experiment
            startTime = experiment[0][1]
            for r in range(0, len(data1)):
                for value in experiment:  # For each [energy-time] value of this experiment
                    if value[0] < ranges1[r]:
                        data1[r].append((value[1] - startTime))  # 0.744
                        break

        print("\nMean: " + str(statistics.median(data1[len(data1)-1])))

        data2 = [[] for x in ranges2]

        for experiment in TimesAndEnergies:  # For each experiment
            startTime = experiment[0][1]
            for r in range(0, len(data2)):
                for value in experiment:  # For each [energy-time] value of this experiment
                    if value[0] < ranges2[r]:
                        if ((value[1] - startTime - 0.01) < 0):
                            data2[r].append(0.01)
                        else:
                            data2[r].append((value[1] - startTime))
                        break

        for D in range(0, 2):
            if D == 0:
                data = copy.deepcopy(data1)
                BoxNames = [str(ranges1[i]) + " to " + str(ranges1[i + 1]) for i in range(0, len(ranges1) - 1)]
                BoxNames.append(" < " + str(ranges1[len(ranges1)-1]))

            else:
                data = copy.deepcopy(data2)
                BoxNames = [str(ranges2[i]) + " to " + str(ranges2[i + 1]) for i in range(0, len(ranges2) - 1)]
                BoxNames.append(" < " + str(ranges2[len(ranges2)-1]))

	    numBoxes = len(BoxNames)
            maxValue = -1
            minValue = 1000
            for box in range(0, numBoxes):
                for d in data[box]:
                    if d < minValue:
                        minValue = d
                    if d > maxValue:
                        maxValue = d

            fig, ax1 = plt.subplots(figsize=(10, 6))
            fig.canvas.set_window_title(name)
            plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

            bp = plt.boxplot(data, notch=0, sym='+', vert=1, whis=1.5)
            plt.setp(bp['boxes'], color='black')
            plt.setp(bp['whiskers'], color='black')
            plt.setp(bp['fliers'], color='red', marker='+')

            # Add a horizontal grid to the plot, but make it very light in color
            # so we can use it for reading data values but not be distracting
            ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
                           alpha=0.5)

            # Hide these grid behind plot objects
            ax1.set_axisbelow(True)
            ax1.set_title('')
            ax1.set_xlabel('Energy (eV)', fontsize=14)
            ax1.set_ylabel('Time (seconds)', fontsize=14)

            # Now fill the boxes with desired colors
            boxColors = ['lightgrey', 'grey', 'red']
            medians = list(range(numBoxes))
            for i in range(numBoxes):
                box = bp['boxes'][i]
                boxX = []
                boxY = []
                for j in range(5):
                    boxX.append(box.get_xdata()[j])
                    boxY.append(box.get_ydata()[j])
                boxCoords = list(zip(boxX, boxY))
                # Alternate between Dark Khaki and Royal Blue
                boxPolygon = Polygon(boxCoords, facecolor=boxColors[0])
                ax1.add_patch(boxPolygon)
                # Now draw the median lines back over what we just filled in
                med = bp['medians'][i]
                medianX = []
                medianY = []
                for j in range(2):
                    medianX.append(med.get_xdata()[j])
                    medianY.append(med.get_ydata()[j])
                    plt.plot(medianX, medianY, 'k')
                    medians[i] = medianY[0]
                # Finally, overplot the sample averages, with horizontal alignment
                # in the center of each box
                plt.plot([np.average(med.get_xdata())], [np.average(data[i])],
                         color='w', marker='', markeredgecolor='k')

            # Set the axes ranges and axes labels
            ax1.set_xlim(0.5, numBoxes + 0.5)
            top = maxValue * 1.5
            bottom = minValue
            if (bottom - 0.01 < 0):
                bottom = 0.01
            bottom = 0.001
            ax1.set_ylim(bottom, top)
            ax1.set_yscale('log')
            ax1.grid(which='minor', color='grey', linestyle='-', linewidth=0.1)
            xtickNames = plt.setp(ax1, xticklabels=np.repeat(BoxNames, 1))

            plt.setp(xtickNames, rotation=45, fontsize=12)

            # Due to the Y-axis scale being different across samples, it can be
            # hard to compare differences in medians across the samples. Add upper
            # X-axis tick labels with the sample medians to aid in comparison
            # (just use two decimal places of precision)

            pos = np.arange(numBoxes) + 1
            upperLabels = [str(np.round(s, 2)) for s in medians]
            weights = ['bold', 'semibold']
            for tick, label in zip(range(numBoxes), ax1.get_xticklabels()):
                k = tick % 2
                ax1.text(pos[tick], top + (top * 0.08), upperLabels[tick],
                         horizontalalignment='center', size='13', weight=weights[k],
                         color=boxColors[2])

            plt.show()




def writeProgressToCSV(progress, procedure_parameters, path, name):
    """
    Stores the execution progress in a csv file.

    progress: list of Progress objects
        Contains all snapshots of the execution.
    procedure_parameters:
        Contains the operations executed. For each operation, there are created five columns in the csv file.
    path: string
        Path to store the csv file
    name: string
        name of the csv file (results_name.csv)
    """

    output_file_text = [[]]
    output_file_text[0] = ["Index"]

    Operations = procedure_parameters[2]

    for i in range(len(Operations)):
        output_file_text[0].append(" ")
        output_file_text[0].append(Operations[i] + " Energy drop (ev/atom)")
        output_file_text[0].append(Operations[i] + " Time (s)")
        output_file_text[0].append(Operations[i] + " Energy (eV/atom)")
        output_file_text[0].append(Operations[i] + " number of loops")
        output_file_text[0].append(Operations[i] + " number of SP Calcs")

    current_operation = -1
    for i in range(len(progress)):
        current_operation = (current_operation + 1) % len(Operations)
        current_phase = progress[i].phase + 1

        if (current_operation == 0):
            output_file_text.append([])
            output_file_text[current_phase].append(str(current_phase))
            output_file_text[current_phase].append(" ")
            output_file_text[current_phase].append("-")
            output_file_text[current_phase].append(str(progress[i].Time))
            output_file_text[current_phase].append(str(progress[i].Energy))
            output_file_text[current_phase].append("-")
            output_file_text[current_phase].append("-")
        else:
            output_file_text[current_phase].append("")
            output_file_text[current_phase].append(str(progress[i].Energy - progress[i-1].Energy))
            output_file_text[current_phase].append(str(progress[i].Time))
            output_file_text[current_phase].append(str(progress[i].Energy))
            output_file_text[current_phase].append(str(progress[i].numberOfLoops))
            output_file_text[current_phase].append(str(progress[i].numberOfSinglePointCalculations))

    csv.register_dialect('myDialect',
                         delimiter=';',
                         quoting=csv.QUOTE_NONE,
                         skipinitialspace=True)
    with open(path + '/results' + name + '.csv', 'w') as f:
        writer = csv.writer(f, dialect='myDialect')
        for row in output_file_text:
            writer.writerow(row)
    f.close()




# =============================================================== #
#                             OPERATIONS                          #
# =============================================================== #


def relaxStructure(s, operation):
    """
    Performs a gulp relaxation (either relaxation of unit cell parameters only, or both unit cell and atomic positions.

    s: structure_class
    operation: string
        Specifies the gulp calculation to execute.

    Returns
    -------
        s : structure_class
            The object containing the relaxed configuration.
        {True or False}:
            Relaxation succeed or failed.
        result: string
            The gulp output.
    """
    s2 = copy.copy(s)
    s.structure, s.energy, result, calc_time = gulp_calc.gulp_relaxation(copy.copy(s), operation, gulp_shells=[])
    return (s, True, result) if result == "converged" else (s2, False, result)




def randomValidStructure(s):
    """
    Computes a random valid configuration of atoms.
    s: structure_class
        The structure_class contains information about the number of atoms, types of atoms and the unit cell parameters.

    Returns
    --------
    s: structure_class
        It contains the valid random configutaion of atoms.
    """
    i = 0
    millis = int(round(time.time() * 1000))
    rng = rngs.NR_Ran(millis)
    rng.warm_up(10)

    say("Trying to find a valid random configuration of atoms...")
    scaled_positions = s.structure.get_scaled_positions()
    print(s.structure.get_cell_lengths_and_angles())
    sys.stdout.write("[%s]" % (" " * len(s.structure.positions)))
    sys.stdout.flush()
    sys.stdout.write("\b" * (len(s.structure.positions) + 1))  # return to start of line, after '['

    while(i<len(s.structure.positions)):
        valid = False
        timer = 10000
        while (not valid):
            scaled_positions[i] = [rng.real_range(0.0, 1.0),
                                   rng.real_range(0.0, 1.0),
                                   rng.real_range(0.0, 1.0)]
            s.structure.set_scaled_positions(scaled_positions)
            valid, overlappingAtoms = s.overlaps(atom=i, check_up_to=i, get_all_overlapping_atoms=False)
            timer -= 1
            if(not valid and timer <= 0):
                i = -1
                say("]\nCould't find a valid random structure. Trying again...\n")
                sys.stdout.write("[")
                break
        i += 1
        sys.stdout.write("-" + str(i))
        sys.stdout.flush()

    sys.stdout.write("]")
    say("\nFound a random valid configuration of atoms.\n")
    return s





# =============================================================== #
#                        ADDITIONAL METHODS                       #
# =============================================================== #

def make_folders(base_path, procedure_parameters):
    """
    Creates all directories needed. For each operation of the optimization function, it creates a folder
    which will contain the .cif files (configurations of atoms).

    base_path: string
        Specifies the path to create the sub-directories.
    procedure_parameters: list
        Contains information about the operations to be executed.

    Returns
    -------
    base_path:
        Returns the path where the sub-directories were created.
    """

    try:
        os.mkdir(base_path)
    except OSError:
        print("Error creating directory " + base_path)
        #sys.exit("Terminating execution: Creation of the directory %s failed" % base_path)

    for operation in procedure_parameters[2]:
        path = ''
        try:
            path = os.path.join(base_path, str(operation))
            os.mkdir(path)
        except OSError:
            sys.exit("Terminating execution: Creation of the directory %s failed" % path)
        else:
            print ("Successfully created the directory %s " % path)
    return base_path


def save_structure(s, path, name, index):
    try:
        s.symmetrize()
        current_structure = cn_plus.Structure_data(s.structure.copy(), 0, s.getEnergy(False), 0.0, [], [])
        current_structure.write_cif(path + '/' + name + str(index) + '.cif')
    except:
        say("Error: Failed to save structure " + path + '/' + name + str(index) + '.cif')


def say(text):
    """
    Prints the input text, and stores it in the output file specified in the .input file.

    text: string
    """
    global output
    print(text)
    output.write(text)


