import numpy as np
import copy
from . import parser
import glob
from os import path

attachment_atom_molecule = 'P'
attachment_atom_crystal = 'R'
crystal_atom = 'Si'
replacement_atom = 'C'



def homogeneous_sphere(h_atoms, percentage, center):
    ''' Generates homogeneous points on a sphere with a radius approximately that of the hydrogens to the center'''
    homogeneous_vectors = []
    samples = int(len(hydrogens) * percentage / 100)
    size = np.linalg.norm([h_atoms[0].x - center.get_x(),h_atoms[0].y - center.get_y(),h_atoms[0].z - center.get_z()])
    phi = np.pi * (3. - np.sqrt(5.))

    # Generate points on a homogeneous sphere
    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2
        theta = phi * i
        x = np.cos(theta)
        z = np.sin(theta)

        vector = parser.Vector(x * size+ center.get_x(), y * size + center.get_y(), z * size + center.get_z())
        homogeneous_vectors.append(vector)

    return homogeneous_vectors

def find_closest_atoms(vectors, atoms):
    ''' Finds the closest atom for every point in a list of points and atoms'''
    attachment_atoms = []
    for point in vectors:
        point_vector = np.array([point.get_x(), point.get_y(), point.get_z()])
        d_min = float('inf')
        for atom in atoms:
            atom_vector = np.array([atom.x,atom.y,atom.z])
            d = np.linalg.norm(atom_vector-point_vector)
            if d < d_min:
                d_min = d
                closest_atom = atom
        closest_atom.kind = 'R'
        attachment_atoms.append(closest_atom)
    return attachment_atoms

if __name__ == '__main__':
    # Generate crystal with the replaced hydrogens
    parser1 = parser.Parser(parser.Reader(2),parser.Interpreter(), parser.Combiner(),  parser.Writer(' '))
    atoms = parser1.reader.read(glob.glob(path.join('input/crystal', '*.xyz'))[0])
    center_vector = parser1.interpreter.find_crystal_center(atoms, 'Si')

    hydrogens = []
    for atom in atoms:
        if atom.kind == 'H':
            hydrogens.append(atom)

    homogeneous_vectors = homogeneous_sphere(hydrogens,20,center_vector)
    find_closest_atoms(homogeneous_vectors, hydrogens)
    parser1.combined = atoms
    parser1.write(path.join('output','crystal','output.xyz'))

    # Combine the crystal with the ligand
    parser2 = parser.Parser(parser.Reader(2),parser.Interpreter(), parser.Combiner(),  parser.Writer(' '))
    parser2.interpret_molecule(glob.glob(path.join('input/molecule','*.xyz'))[0],attachment_atom_molecule)
    parser2.interpret_crystal(glob.glob(path.join('output','crystal', '*.xyz'))[0], attachment_atom_crystal,crystal_atom,replacement_atom)
    parser2.combine()
    parser2.write(path.join('output','output_percentage_coverage.xyz'))