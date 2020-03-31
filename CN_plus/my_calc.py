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
import threading
import gulp_calc
import time
import inputs
import itertools as it
import multiprocessing as mp
from functools import partial
import functions
import logging

# ============================================================================ #
#                                  MY CALC                                     #
# ============================================================================ #


logging.basicConfig(level=logging.DEBUG,
                    format='(%(threadName)-10s) %(message)s',
                    )


class EnergyCalc(object):
    def __init__(self):
        self.lock = threading.Lock()
        self.Coulomb = 0
        self.Buck = 0

    def increment(self, coulomb, buck):
        self.lock.acquire()
        try:
            self.Coulomb = self.Coulomb + coulomb
            self.Buck = self.Buck + buck
        finally:
            self.lock.release()


class matrix:
    depth = 2
    DipoleConstant = 30.15080016
    C = 14.399645351950543
    atoms_data = ""

    def __init__(self, structure, atoms_data, numberOfAtoms, depth=2):
        self.buckingham = [[0 for i in range(numberOfAtoms)] for j in range(numberOfAtoms)]
        self.atoms_data = atoms_data
        self.charges = [0 for i in range(numberOfAtoms)]
        self.depth = depth
        self.atomTypes = []
        self.symbols = []
        self.charges = []
        self.numberOfAtomsPerType = []
        for entry in atoms_data:
            self.atomTypes.append(entry.split(' ')[0])
            self.numberOfAtomsPerType.append(int(entry.split(' ')[1]))
            for i in range(0, self.numberOfAtomsPerType[len(self.numberOfAtomsPerType) - 1]):
                self.charges.append(int(entry.split(' ')[2]))
                # self.symbols.append(entry.split(' ')[0])

        self.structure = structure
        self.numberOfAtoms = numberOfAtoms

        self.vecs = None  # np.array(direction, triple), written by rows
        self.ions = None  # np.array()
        self.N = None

        self.ions = []
        for pos in self.structure.get_scaled_positions():
            self.ions.append(list(np.concatenate(([0], pos))))

        for i in range(0, len(self.charges)):
            self.ions[i][0] = self.charges[i]

        self.ions = np.array(self.ions)
        self.N = len(self.ions)
        self.vecs = self.structure.get_cell()

    def getEnergy(self):
        return (self.Coulomb(self.depth) - self.DipoleConstant * np.linalg.norm(
            self.dipole_moment_in_A()) ** 2 / self.cell_volume()) / self.numberOfAtoms

    def dipole_moment(self):

        moment = np.zeros(3)
        for ion in self.ions:
            moment += ion[0] * ion[1:]

        return moment

    def dipole_moment_in_A(self):
        return np.matmul(self.dipole_moment(), self.vecs)

    def cell_volume(self):
        return abs(np.linalg.det(self.vecs))

    def Coulomb(self, depth=1):

        runParallel = False

        if (runParallel):
            '''
            # SOLUTION 1 - SOMETHING IS WRONG...
            pool = mp.Pool(mp.cpu_count())
            resultsCoulomb = [pool.apply(functions.myCalc, args=(i, j, depth, self.vecs)) for i, j in it.combinations(self.ions, 2)]
            print(resultsCoulomb)
            energy = sum(resultsCoulomb)
            total_buckingham_energy = 0 # sum(resultsBuckingham)
            pool.close()

            print(energy)
            energy += 0.5 * sum([i[0] ** 2 for i in self.ions]) * self.Coulomb_self(depth)
            print(energy)
            '''

            '''
            # SOLUTION 2 - WORKS (SLOW)
            outputCoulomb = mp.Queue()
            outputBuckingham = mp.Queue()
            processes = [mp.Process(target=self.myCalc, args=(i, j, depth, outputCoulomb, outputBuckingham)) for i, j in it.combinations(self.ions, 2)]
            for p in processes:
                p.start()

            for p in processes:
                p.join()

            resultsCoulomb = [outputCoulomb.get() for p in processes]
            resultsBuckingham = [outputBuckingham.get() for p in processes]

            energy = sum(resultsCoulomb)
            total_buckingham_energy = sum(resultsBuckingham)

            energy += 0.5 * sum([i[0] ** 2 for i in self.ions]) * self.Coulomb_self(depth)
            # print(f'With repulsion {energy}')
            '''

            # SOLUTION 3 - Threads

            # outputCoulomb = mp.Queue()
            # outputBuckingham = mp.Queue()

            energyClass = EnergyCalc()

            temp1 = time.time()

            for i, j in it.combinations(self.ions, 2):
                threading.Thread(target=self.myCalc, args=(i, j, depth, energyClass)).start()

            main_thread = threading.currentThread()

            for t in threading.enumerate():
                if t is not main_thread:
                    t.join()

            print("Time threads: " + str(time.time() - temp1))

            # resultsCoulomb = [outputCoulomb.get() for t in threads]
            # resultsBuckingham = [outputBuckingham.get() for t in threads]

            energy = energyClass.Coulomb
            total_buckingham_energy = energyClass.Buck

            temp1 = time.time()
            energy += 0.5 * sum([i[0] ** 2 for i in self.ions]) * self.Coulomb_self(depth)
            print("Time Coulomb.self: " + str(time.time() - temp1))


        else:
            energy = 0
            total_buckingham_energy = 0
            for i, j in it.combinations(self.ions, 2):
                coulomb, buckingham = self.Coulomb_and_Buckingham(i[1:], j[1:], i[0], j[0], depth)
                energy += i[0] * j[0] * coulomb
                total_buckingham_energy += buckingham
            energy += 0.5 * sum([i[0] ** 2 for i in self.ions]) * self.Coulomb_self(depth)

        return energy + total_buckingham_energy

    def myCalc(self, i, j, depth, energyClass):
        # logging.debug('Started')
        coulomb, buckingham = self.Coulomb_and_Buckingham(i[1:], j[1:], i[0], j[0], depth)
        energyClass.increment(i[0] * j[0] * coulomb, buckingham)
        # outputCoulomb.put(energy)
        # outputBuckingham.put(total_buckingham_energy)
        # return energy
        # logging.debug('Finished')

    def Coulomb_and_Buckingham(self, t1, t2, c1, c2, depth=2):
        p1 = np.matmul(t1, self.vecs)
        p2 = np.matmul(t2, self.vecs)

        coulomb_energy = 0
        buckingham_energy = 0
        specie1 = ""
        specie2 = ""

        if c1 == -2:
            specie1 = "O"
            rad1 = 1.21
        elif c1 == 2:
            specie1 = "Sr"
            rad1 = 1.32
        else:
            specie1 = "Ti"
            rad1 = 0.56

        if c2 == -2:
            specie2 = "O"
            rad2 = 1.21
        elif c2 == 2:
            specie2 = "Sr"
            rad2 = 1.32
        else:
            specie2 = "Ti"
            rad2 = 0.56

        constantA, constantB, constantC = self.getBuck(specie1, specie2)

        for shift in it.product(range(depth, -depth - 1, -1), repeat=3):
            b = np.matmul(np.array(shift), self.vecs)
            # b = self.myMul(shift, self.vecs) # np.matmul(np.array(shift), self.vecs)
            r = self.myNorm(p2 + b - p1)
            # r = np.linalg.norm(b) # Too slow...
            coulomb_energy += self.C / r
            if constantB != 0:
                buckingham_energy += constantA * math.exp(-r / constantB) - constantC / (r ** 6)
            else:
                buckingham_energy += constantC / (r ** 6)

        return coulomb_energy, buckingham_energy

    def myMul(self, a, b):
        return [a[0] * b[0][0] + a[1] * b[0][1] + a[2] * b[0][2], a[0] * b[1][0] + a[1] * b[1][1] + a[2] * b[1][2],
                a[0] * b[2][0] + a[1] * b[2][1] + a[2] * b[2][2]]

    def myNorm(self, l):
        s = 0
        for i in l:
            s += i ** 2
        return math.sqrt(s)

    def Coulomb_self(self, depth=2):
        energy = 0
        C = 14.399645351950543
        shifts = list(it.product(range(depth, -depth - 1, -1), repeat=3))
        shifts.remove((0, 0, 0))

        for shift in shifts:
            r = np.linalg.norm(np.matmul(np.array(shift), self.vecs))
            energy += C / r

        return energy

    def getBuck(self, atom1="O", atom2="O"):
        if atom1 == "O" and atom2 == "O":
            return 1388.77, 0.36262, 175
        elif (atom1 == "O" and atom2 == "Sr") or (atom1 == "Sr" and atom2 == "O"):
            return 1952.390, 0.33685, 19.22000
        elif (atom1 == "O" and atom2 == "Ti") or (atom1 == "Ti" and atom2 == "O"):
            return 4590.7279, 0.261000, 0.0000000
        else:
            return 0, 0, 0

    def Ewald(self, alpha=-1):
        '''
        Large alpha forces real part to converge faster. Normally, it is around 5/L
        '''
        if alpha < 0:
            alpha = 2 / self.vecs[0, 0]
        realDepth = 4
        reciprocalDepth = 5
        energy = 0
        # np.seterr(all='warn', over='raise')

        # Real part
        # for i,j in it.combinations(self.ions, 2):
        #     # print(self.__ewald_real_coeff(alpha, i[1:], j[1:], 3))
        #     energy += i[0]*j[0]*self.__ewald_real_coeff(alpha, i[1:], j[1:], 3)

        for i_index, i in enumerate(self.ions):
            for j_index, j in enumerate(self.ions):
                if i_index != j_index:
                    r = np.linalg.norm(np.matmul(i[1:], self.vecs) - np.matmul(j[1:], self.vecs))
                    energy += i[0] * j[0] * math.erfc(alpha * r) / r

                shifts = list(it.product(range(realDepth, -realDepth - 1, -1), repeat=3))
                shifts.remove((0, 0, 0))

                for shift in shifts:
                    r = np.linalg.norm(
                        np.matmul(i[1:], self.vecs) + np.matmul(np.array(shift), self.vecs) - np.matmul(j[1:],
                                                                                                        self.vecs))
                    energy += i[0] * j[0] * math.erfc(alpha * r) / r

        energy = energy / 2

        print("Ewald, Real part:", energy)
        # self interaction term
        energy_self = -alpha * sum([i[0] ** 2 for i in self.ions]) / math.sqrt(math.pi)
        # print("Ewald, Self interaction term:", energy_self)
        energy += energy_self

        # Ewald reciprocal space
        energy_reciprocal = 0

        for i_index, i in enumerate(self.ions):
            for j_index, j in enumerate(self.ions):

                shifts = list(it.product(range(reciprocalDepth, -reciprocalDepth - 1, -1), repeat=3))
                shifts.remove((0, 0, 0))

                for shift in shifts:
                    # NOTE, here k is defined assuming cubic cell, I have no idea how k should be defined for non-cubic cells
                    k = np.array(shift) * 2 * math.pi / self.vecs[0, 0]
                    term = i[0] * j[0]
                    term = (term * 4 * math.pi ** 2) / np.matmul(k, k)
                    term = term * math.exp(-np.matmul(k, k) / (4 * alpha ** 2))
                    r = np.matmul(j[1:], self.vecs) - np.matmul(i[1:], self.vecs)
                    # print(k, r)
                    term = term * math.cos(np.matmul(k, r))
                    energy_reciprocal += term

        energy_reciprocal = energy_reciprocal / (2 * math.pi * self.cell_volume())
        print("Ewald, reciprocal part:", energy_reciprocal + energy_self)

        energy += energy_reciprocal
        # Unit conversion
        energy = energy * self.C
        return energy / self.numberOfAtoms

    def updateMatrix(self):
        time1 = time.time()
        print("Matrix for depth: " + str(self.depth))

        tmp_structure = self.structure.copy()
        scaled_positions = tmp_structure.get_scaled_positions().copy()
        timeA = time.time()
        timeB = 0

        for i in range(0, self.numberOfAtoms):
            for j in range(0, self.numberOfAtoms):
                # print("(i,j) = " + str(i) + ", " + str(j))
                self.coulomb[i][j] = 0
                self.buckingham[i][j] = 0
                pos1 = self.structure.get_scaled_positions()[i]
                pos2 = self.structure.get_scaled_positions()[j]
                scaled_positions[0] = copy.copy(pos1)

                for x in range(-self.depth, self.depth + 1):  # [-depth, 0, depth]
                    for y in range(-self.depth, self.depth + 1):
                        for z in range(-self.depth, self.depth + 1):
                            scaled_positions[1] = [pos2[0] + x, pos2[1] + y, pos2[2] + z]
                            tmp_structure.set_scaled_positions(scaled_positions)
                            timetmp = time.time()
                            distance = tmp_structure.get_distance(0, 1)
                            # print("dist " + str(i) + ", " + str(j) + " - (" + str(x) + ", " + str(y) + "," + str(z) + "): " + str(distance))
                            timeB += time.time() - timetmp
                            if (distance > 0):
                                if (i == j and x == 0 and y == 0 and z == 0):
                                    continue
                                else:
                                    self.coulomb[i][j] += self.C / distance
                                    A, B, C = self.getBuck(self.symbols[i], self.symbols[j])
                                    if (B != 0):
                                        self.buckingham[i][j] += A * math.exp(-distance / B) - C / (distance ** 6)
        # print("Total time: " + str(time.time() - timeA))
        # print("Time get distance: " + str(timeB))

    def getEnergy2(self):

        self.coulomb = [[0 for i in range(self.numberOfAtoms)] for j in range(self.numberOfAtoms)]
        self.buckingham = [[0 for i in range(self.numberOfAtoms)] for j in range(self.numberOfAtoms)]
        self.charges = [0 for i in range(self.numberOfAtoms)]
        self.depth = self.depth
        self.charges = []
        self.numberOfAtomsPerType = []
        for entry in self.atoms_data:
            self.atomTypes.append(entry.split(' ')[0])
            self.numberOfAtomsPerType.append(int(entry.split(' ')[1]))
            for i in range(0, self.numberOfAtomsPerType[len(self.numberOfAtomsPerType) - 1]):
                self.charges.append(int(entry.split(' ')[2]))
                self.symbols.append(entry.split(' ')[0])

        self.updateMatrix()

        coulomb_m = np.matrix(self.coulomb)
        charges_m = np.matrix(self.charges)
        tmp1 = np.matmul(charges_m, coulomb_m)
        Coulomb = np.matmul(tmp1, np.transpose(charges_m))

        Buckingham = 0
        for i in range(0, self.numberOfAtoms):
            for j in range(0, self.numberOfAtoms):
                Buckingham += self.buckingham[i][j]

        DipoleConstant = 30.15080016

        self.ions = []
        for pos in self.structure.get_scaled_positions():
            self.ions.append(list(np.concatenate(([0], pos))))

        for i in range(0, len(self.charges)):
            self.ions[i][0] = self.charges[i]

        self.ions = np.array(self.ions)

        dipoleMoment = np.zeros(3)
        for ion in self.ions:
            dipoleMoment += ion[0] * ion[1:]

        CoulombCorrected = 0.5 * Coulomb.item(0) - self.DipoleConstant * np.linalg.norm(
            self.dipole_moment_in_A()) ** 2 / self.cell_volume()
        totalEnergy = (CoulombCorrected + Buckingham) / self.numberOfAtoms
        # print("Coulomb corrected: " + str(CoulombCorrected/self.numberOfAtoms))
        # print("Buck: " + str(Buckingham/self.numberOfAtoms))
        # print("Coulomb + buck: " + str((CoulombCorrected + Buckingham)/self.numberOfAtoms))
        # print("With constant: " + str((constant*CoulombCorrected + Buckingham)/self.numberOfAtoms))
        return totalEnergy

