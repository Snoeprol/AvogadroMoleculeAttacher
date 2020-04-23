import numpy as np
import copy
import parser as ps
import glob
from os import path
def check_powerset_entities(powerset_list):
    """ Checks how many 1's a subset in the form [[0], [1], [0]] etc. contains """
    entities = 0
    for i in powerset_list:
        if i == 1:
            entities += 1
    return entities
    
def generate_powerset_subsets(size, subset_size):
    iteration_list = []
    sets = []
    for i in range(size):
        iteration_list.append(0)
    if check_powerset_entities(iteration_list) == subset_size:
        sets.append(copy.deepcopy(iteration_list))
    
    # Generate subsets of length subset_size
    for i in range(2**size):
        for l in range(size):
            if iteration_list[l] + 1 == 1:
                iteration_list[l] += 1
                for k in range(l):
                    iteration_list[k] = 0
                if check_powerset_entities(iteration_list) == subset_size:
                    sets.append(copy.deepcopy(iteration_list))
                break
            else:
                continue
    return sets

def calc_average(positions):
    total = 0
    for position in positions:
        total += position
    return total / len(positions)

def calculate_quasi_st_dev(atoms):
    """ Spits out a number that is a measurement of the spread """
    x_positions = [atom.x for atom in atoms]
    y_positions = [atom.y for atom in atoms]
    z_positions = [atom.z for atom in atoms]

    average_x = calc_average(x_positions)
    average_y = calc_average(y_positions)
    average_z = calc_average(z_positions)

    sum = 0
    for atom in atoms:
        sum += (atom.x - average_x)**2 + (atom.y - average_y)**2 +(atom.z - average_z)**2
    return sum

def find_homogeneous_points(atoms, percentage):

    atom_count = len(atoms)
    amount_in_subset = int(float(percentage) * atom_count/ 100)
    quasi_st_devs = []
    subsets = generate_powerset_subsets(atom_count,amount_in_subset)

    # Calculate quasi-st-devs
    for subset in subsets:
        atoms_subset = []
        for atom_index, bit in enumerate(subset):
            if bit == 1:
                atoms_subset.append(atoms[atom_index])
        
        quasi_st_devs.append(calculate_quasi_st_dev(atoms_subset))
    
    best_subset_index = max(range(len(quasi_st_devs)), key=quasi_st_devs.__getitem__)
    best_subset_atoms = []
    for index, bit in enumerate(subsets[best_subset_index]):
        if bit == 1:
            best_subset_atoms.append(atoms[index])
    return best_subset_atoms


parser1 = ps.Parser(ps.Reader(2),ps.Interpreter(), ps.Combiner(),  ps.Writer(' '))
atoms = parser1.reader.read(glob.glob(path.join('input/molecule', '*.xyz'))[0])
print(find_homogeneous_points(atoms[0:10], 20))