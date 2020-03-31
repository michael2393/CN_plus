try:
    from ase.calculators.gulp import GULP
except ImportError:
    from ase.calculators.gulp import Gulp as GULP

import spglib as sp
import ase

# ============================================================================ #
#                                 STRUCTURE                                    #
# ============================================================================ #

class structure_class:
    """
    This class contains all the information needed in an 'optimization' procedure.
    It stores the configuration of atoms, the types, number, charge, and radius of them,
    and also stores infomation such as the number of single point calculations performed.
    Finally, it contains all necessary methods to handle the configuration.
    """

    def __init__(self, params):
        self.params = params
        self.grid_points = params["grid_points"]["value"]
        self.relax_structure_keywords = params["gulp_keywords_relax"]["value"]
        self.relax_structure_options = params["gulp_options_relax"]["value"]

        self.relax_unit_cell_keywords = params["gulp_keywords_relax_unit_cell"]["value"]
        self.relax_unit_cell_options = params["gulp_options_relax_unit_cell"]["value"]

        #self.relaxation_gnorms = params["gulp_max_gnorm"]["value"]
        self.max_gnorm = params["gulp_max_gnorm"]["value"]
        self.gulp_library = params["gulp_library"]["value"]
        print(self.gulp_library)
        self.atoms_data = []
        self.distinctAtoms = []     # [symbol, number, charge, radius]
        self.Atoms = []             # [symbol, charge, radius]
        self.structure = []
        self.energy = 0
        self.numberOfAtoms = 0
        self.numberOfSinglePointCalculations = 0


    def initializeAtomsData(self):
        print(self.atoms_data)
        for species_data in self.atoms_data:
            # Extract species, number and charge of atoms
            specie, num_atoms, ionic_charge, radius = species_data.split()
            num_atoms = int(num_atoms)
            self.numberOfAtoms += num_atoms
            self.distinctAtoms.append({ "symbol": specie,
                                        "numberOfAtoms": int(num_atoms),
                                        "charge": float(ionic_charge),
                                        "radius": float(radius)})
            for i in range(0, int(num_atoms)):
                self.Atoms.append({ "symbol": specie,
                                    "charge": float(ionic_charge),
                                    "radius": float(radius)})

    def getCharges(self):
        species, num_atoms, ionic_charge, radius = self.atoms_data.split()
        return ionic_charge

    def getNumberOfAtoms(self):
        species, num_atoms, ionic_charge, radius = self.atoms_data.split()
        return ionic_charge

    def getSpecies(self):
        species, num_atoms, ionic_charge, radius = self.atoms_data.split()
        return species

    def getRadii(self):
        species, num_atoms, ionic_charge, radius = self.atoms_data.split()
        return ionic_charge

    def getEnergy(self, recalculate):
        """
        It calculates the energy of the unit cell.

        recalculate: boolean
            If True, it performs a new single point calculation (using gulp) to re-calculate the energy.
            If False, it return the last known energy calculated in a previous point.

        Return
        ------
        It returns a boolean which specifies if the returned energy is correct or not (due to gulp failure, and
        the average energy per atom.
        """

        if recalculate:
            try:
                calc = (GULP(keywords=self.params["gulp_keywords_calculate_energy"]["value"],
                             library=self.params["gulp_library"]["value"]))
                self.structure.set_calculator(calc)
                self.energy = self.structure.get_potential_energy()
                self.numberOfSinglePointCalculations += 1
                return True, self.energy / self.numberOfAtoms
            except:
                return False, 0

        else:
            return True, self.energy / self.numberOfAtoms


    def overlaps(self, atom, check_up_to, get_all_overlapping_atoms=True):
        """
        This method checks whether a specific atom overlaps with any other atoms in the unit cell.

        atom: int
            The atom to check whether it overlaps with other atoms.
        check_up_to: int
            This method checks if the given atom overlaps with the rest of the atoms with indexes 0 to 'check_up_to'.
        get_all_overlapping_atoms: boolean
            If this parameter is True, then this method searches for all possible overlaps, otherwise, if
            one overlapping atoms is found, it stops.

        Return:
        ------
        valid: boolean
            This parameter specifies if the given atom is 'valid' or not
        overlappingAtoms: list
            It contains all indexes with the atoms that the given atom overlaps with.
        """
        if (check_up_to == 0):
            return True, []
        distances = self.structure.get_distances(atom, [i for i in range(0, check_up_to)], mic=True)
        minimum_percentage_allowed = 0.99
        valid = True
        overlappingAtoms = []

        init_distance = self.Atoms[atom]["radius"]

        for i in range(0, check_up_to):
            if (i == atom):
                continue
            minimum_distance = init_distance + self.Atoms[i]["radius"]
            if (distances[i] / minimum_distance < minimum_percentage_allowed):
                overlappingAtoms.append([i, minimum_distance - distances[i]])
                #print("Minimum allowed: " + str(minimum_distance) + ", dist: " + str(distances[i]))
                valid = False
                if (not get_all_overlapping_atoms):
                    break

        return valid, overlappingAtoms


    def symmetrize(self):
        try:
            lattice, scaled_pos, atomic_numbers = sp.standardize_cell(self.structure)
            self.structure = ase.Atoms(cell=lattice, scaled_positions=scaled_pos, numbers=atomic_numbers,
                                          charges=self.structure.get_initial_charges(),
                                          pbc=[True, True, True])
            print("Symmetrize atoms: success")
        except:
            print("Symmetrize atoms: error")
