import numpy as np
import copy
from . import parser
import glob
from os import path
import tkinter as tk

def homogeneous_sphere(h_atoms, percentage, center):
    ''' Generates homogeneous points on a sphere with a radius approximately that of the hydrogens to the center'''
    homogeneous_vectors = []
    samples = int(len(h_atoms) * percentage / 100)
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

def main(attachment_atom_molecule = 'P', crystal_atom = 'Si',replacement_atom = 'C',percentage = 10):
    '''Replaces percentage of hydrogens on input crystal and replaces it with input molecule'''

    # Generate crystal with the replaced hydrogens
    parser1 = parser.Parser(parser.Reader(2),parser.Interpreter(), parser.Combiner(),  parser.Writer(' '))
    atoms = parser1.reader.read(glob.glob(path.join('input/crystal', '*.xyz'))[0])
    center_vector = parser1.interpreter.find_crystal_center(atoms, 'Si')

    hydrogens = []
    for atom in atoms:
        if atom.kind == 'H':
            hydrogens.append(atom)
    homogeneous_vectors = homogeneous_sphere(hydrogens,percentage,center_vector)
    find_closest_atoms(homogeneous_vectors, hydrogens)
    parser1.combined = atoms
    parser1.write(path.join('output','output.xyz'))

    # Combine the crystal with the ligand
    attachment_atom_crystal = 'R'
    parser2 = parser.Parser(parser.Reader(2),parser.Interpreter(), parser.Combiner(),  parser.Writer(' '))
    parser2.interpret_molecule(glob.glob(path.join('input/molecule','*.xyz'))[0],attachment_atom_molecule)
    parser2.interpret_crystal(glob.glob(path.join('output','crystal', '*.xyz'))[0], attachment_atom_crystal,crystal_atom,replacement_atom)
    parser2.combine()
    parser2.write(path.join('output','output.xyz'))

def generate():
    main(e1.get(), e2.get(), e3.get(), int(e4.get()))

if __name__ == '__main__':
    ''' Interface '''
    master = tk.Tk()
    master.title("Meerio Crystal Attacher")
    tk.Label(master, 
            text="Molecule Attachment Atom (e.g. : P):").grid(row=0)
    tk.Label(master, 
            text="Crystal Atom (e.g. : Si):").grid(row=1)
    tk.Label(master,
            text = "Replacement Atom for Attachment atom (e.g. : O):").grid(row=2)
    tk.Label(master,
            text = "Percentage of Hydrogen that should be replaced (e.g. : 10)").grid(row=3)

    e1 = tk.Entry(master)
    e2 = tk.Entry(master)
    e3 = tk.Entry(master)
    e4 = tk.Entry(master)
    e1.grid(row=0, column=1)
    e2.grid(row=1, column=1)
    e3.grid(row=2, column=1)
    e4.grid(row=3, column=1)

    tk.Button(master, 
            text='Quit', 
            command=master.quit).grid(row=5, 
                                        column=0, 
                                        sticky=tk.W, 
                                        pady=4)
    tk.Button(master, 
            text='Generate', command=generate).grid(row=5, 
                                                        column=1, 
                                                        sticky=tk.W, 
                                                        pady=4)
    tk.mainloop()
