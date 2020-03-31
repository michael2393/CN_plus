try:
    from ase.calculators.gulp import GULP
except ImportError:
    from ase.calculators.gulp import Gulp as GULP

import numpy as np
import gulp_calc
import time as T
import random
import copy
import itertools
import rngs
import time

# =============================================================================
# =============================================================================
def axisHeuristic(s):
    """
    This function implements a greedy strategy and tries to minimize the energy by searching the 'axes' neighbourhood.

    s: structure_class
        Contains the configuration of atoms.

    Returns
    -------
    s: structure_class
    boolean:
        True if managed to find a better (lower energy) structure, otherwise False
    numberOfLoops:
        The maximum number of times that an atom changed position.
    numberOfSinglePointCalculations: int
        Number of energy calculations.
    Energies_Times: list of lists [[energy_1, time_1], [energy_2, time_2], ...]
        This list contains the times that the local search reached the corresponding energies.
    """

    grid_points = s.grid_points
    steps = [s.structure.get_cell_lengths_and_angles()[0] / grid_points[0],
             s.structure.get_cell_lengths_and_angles()[1] / grid_points[1],
             s.structure.get_cell_lengths_and_angles()[2] / grid_points[2]]


    improved = True
    iteration = 0
    numberOfLoops = 0
    timeout = s.numberOfAtoms
    initial_energy = s.getEnergy(True)
    print("Initial energy: " + str(s.getEnergy(True)))

    Energies_Times = []

    while improved:
        atom = iteration % s.numberOfAtoms
        numberOfLoops += 1 if (atom == 0) else 0

        initialPosition = copy.copy(s.structure[atom].position)
        bestPosition = copy.copy(s.structure[atom].position)
        success, tmp_energy = s.getEnergy(True)
        Energies_Times.append([tmp_energy, time.time()])
        best_energy = copy.copy(tmp_energy)

        overlapsResolved = False


        for axes in range(0,3):
            s.structure[atom].position = copy.copy(initialPosition)

            for i in range(0, grid_points[axes]):
                s.structure[atom].position[axes] = i * steps[axes]
                valid, overlappingAtoms = s.overlaps(atom, s.numberOfAtoms, get_all_overlapping_atoms=False)

                if(valid):
                    success, current_energy = s.getEnergy(True)
                    if(current_energy < best_energy - 0.0):
                        bestPosition = copy.copy(s.structure[atom].position)
                        best_energy = current_energy
                        Energies_Times.append([best_energy, time.time()])
                        print("Improved: " + str(best_energy) + ", resolve overlaps: " + str(overlapsResolved))
                        timeout = s.numberOfAtoms
                        overlapsResolved = False

        if(best_energy >= tmp_energy):
            s.structure[atom].position = initialPosition
        else:
            s.structure[atom].position = bestPosition


        iteration += 1
        timeout -= 1
        improved = False if timeout <= 0 else True


    print("Energy after axes: " + str(s.getEnergy(True)))

    return s, True if (s.getEnergy(True) < initial_energy) else False, numberOfLoops, s.numberOfSinglePointCalculations, Energies_Times





def swapAllHeuristic(s):
    """
    This function implements a greedy strategy and tries to minimize the energy by searching the 'swap all' neighbourhood'.

    s: structure_class
        Contains the configuration of atoms.

    Returns
    -------
    s: structure_class
    boolean:
        True if managed to find a better (lower energy) structure, otherwise False
    numberOfLoops:
        The maximum number of times that an atom changed position.
    numberOfSinglePointCalculations: int
        Number of energy calculations.
    """

    grid_points = s.grid_points
    steps = [s.structure.get_cell_lengths_and_angles()[0] / grid_points[0],
             s.structure.get_cell_lengths_and_angles()[1] / grid_points[1],
             s.structure.get_cell_lengths_and_angles()[2] / grid_points[2]]

    positionsEachAxis = [[i * steps[0] for i in range(0, grid_points[0] + 1)],
                 [i * steps[1] for i in range(0, grid_points[1] + 1)],
                 [i * steps[2] for i in range(0, grid_points[2] + 1)]]

    positions = list(itertools.product(*positionsEachAxis))

    improved = True
    iteration = 0
    numberOfLoops = 0
    timeout = s.numberOfAtoms
    initial_energy = s.getEnergy(True)
    print("Initial energy: " + str(s.getEnergy(True)))

    while improved:
        atom = iteration % s.numberOfAtoms
        numberOfLoops += 1 if (atom == 0) else 0

        initialPosition = copy.copy(s.structure[atom].position)
        bestPosition = copy.copy(s.structure[atom].position)
        success, tmp_energy = s.getEnergy(True)
        best_energy = copy.copy(tmp_energy)

        overlapsResolved = False

        for pos in positions:
            s.structure[atom].position = pos
            valid, overlappingAtoms = s.overlaps(atom, s.numberOfAtoms,
                                                 get_all_overlapping_atoms=False)  # True if we want to resolve them

            if(valid):
                success, current_energy = s.getEnergy(True)
                if(current_energy < best_energy - 0.0):
                    bestPosition = copy.copy(s.structure[atom].position)
                    best_energy = current_energy
                    print("Improved: " + str(best_energy) + ", resolve overlaps: " + str(overlapsResolved))
                    timeout = s.numberOfAtoms
                    overlapsResolved = False

        if(best_energy >= tmp_energy):
            s.structure[atom].position = initialPosition
        else:
            s.structure[atom].position = bestPosition


        iteration += 1
        timeout -= 1
        improved = False if timeout <= 0 else True


    print("Energy after swap all pairs heuristic: " + str(s.getEnergy(True)))

    return s, True if (s.getEnergy(True) < initial_energy) else False, numberOfLoops, s.numberOfSinglePointCalculations




def swapAtomsHeuristic(s):
    """
    This function implements a greedy strategy and tries to minimize the energy by searching the 'swap atoms' neighbourhood'.

    s: structure_class
        Contains the configuration of atoms.

    Returns
    -------
    s: structure_class
    boolean:
        True if managed to find a better (lower energy) structure, otherwise False
    numberOfLoops:
        The maximum number of times that an atom changed position.
    numberOfSinglePointCalculations: int
        Number of energy calculations.
    """

    improved = True
    iteration = 0
    numberOfLoops = 0
    timeout = s.numberOfAtoms
    initial_energy = s.getEnergy(True)
    print("Initial energy: " + str(s.getEnergy(True)))


    while improved:
        atom = iteration % s.numberOfAtoms
        numberOfLoops += 1 if (atom == 0) else 0

        success, tmp_energy = s.getEnergy(True)
        best_energy = copy.copy(tmp_energy)

        overlapsResolved = False

        for i in range(0, s.numberOfAtoms):
            if(s.Atoms[atom]["symbol"] == s.Atoms[i]["symbol"]):
                continue
            s = swap(s, atom, i)
            valid, overlappingAtoms = s.overlaps(atom, s.numberOfAtoms,get_all_overlapping_atoms=False)  # True if we want to resolve them
            valid2, overlappingAtoms = s.overlaps(i, s.numberOfAtoms, get_all_overlapping_atoms=False) # True if we want to resolve them

            if(valid and valid2):
                success, current_energy = s.getEnergy(True)
                if(current_energy < best_energy - 0.0):
                    best_energy = current_energy
                    print("Improved: " + str(best_energy) + ", resolve overlaps: " + str(overlapsResolved))
                    timeout = s.numberOfAtoms
                    overlapsResolved = False
                else:
                    s = swap(s, atom, i)
            else:
                s = swap(s, atom, i)

        iteration += 1
        timeout -= 1
        improved = False if timeout <= 0 else True


    print("Energy after swap atoms: " + str(s.getEnergy(True)))

    return s, True if (s.getEnergy(True) < initial_energy) else False, numberOfLoops, s.numberOfSinglePointCalculations



def swap(s, i, j):
    """
    Swaps the positions of the given atoms.

    s: structure_class
    i: int
        The index of the first atom.
    j: int
        The index of the second atom.

    Returns
    -------
    s: structure_class
    """

    tmpPos = copy.copy(s.structure[j].position)
    s.structure[j].position = copy.copy(s.structure[i].position)
    s.structure[i].position = copy.copy(tmpPos)
    return s

