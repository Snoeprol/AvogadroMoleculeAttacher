from os import path
import numpy as np
import glob
from . import mathhelper as mh
import copy

class Atom:

    def __init__(self, kind, x, y, z):
        self.x = x 
        self.y = y
        self.z = z           
        self.kind = kind

class Vector:

    def __init__(self, x, y, z):
        self.vector = np.array([x, y,z])

    def get_x(self):
        return self.vector[0]

    def get_y(self):
        return self.vector[1]
    
    def get_z(self):
        return self.vector[2]

    def normalize(self):
        norm = (self.get_x()**2 + self.get_y()**2 + self.get_z()**2)**0.5
        self.vector[0] /= norm
        self.vector[1] /= norm
        self.vector[2] /= norm

class Molecule:

    def __init__(self, atoms):
        self.atoms = atoms

class Parser: 

    def __init__(self, reader, interpreter, combiner, writer):
        self.reader = reader
        self.interpreter = interpreter
        self.combiner = combiner
        self.writer = writer
        self.output = []
        
    def interpret_crystal(self, input_path, attachment_atom, crystal_atom, replacement_atom):
        """ Interpreted content is a list of atoms, attachment atoms, center and normal vectors"""
        content = self.reader.read(input_path)
        interpreted_content = self.interpreter.interpret(content, attachment_atom, crystal_atom, replacement_atom)
        self.output.append(interpreted_content)
        
    def interpret_molecule(self,input_path, attachment_atom):
        """ Interpreted content is a list of atoms and the center of mass vector"""
        content = self.reader.read(input_path)
        interpreted_content = self.interpreter.interpret(content, attachment_atom)
        self.output.append(interpreted_content)

    def combine(self):
        """Combines the data"""
        self.combined = self.combiner.combine(self.output)

    def write(self, output_path):
        """Writes the data to output_path"""
        self.writer.write(self.combined, output_path)
    
class Reader:

    def __init__(self,lines_to_skip = 0):
        self.lines_to_skip = lines_to_skip

    def read(self,input_path):
        content = self.read_file(input_path)
        if '.xyz' in input_path:
            return self.map_lines_to_atoms(content)
        elif '.csv' in input_path:
            return self.map_lines_to_elements(content)
        
    def read_file(self, input_path):
        with open(input_path, 'r') as stream:
            return stream.readlines()[self.lines_to_skip:]

    def map_lines_to_atoms(self,lines):
        atoms = []
        for line in lines:
            atom = self.map_line_to_atom(line)
            if atom is not None:
                atoms.append(self.map_line_to_atom(line))
        return atoms

    def map_line_to_atom(self,line):
        cleaned_line = self.clean_line(line)
        if cleaned_line[0] == '\n':
            return
        else:
            return Atom(cleaned_line[0] , float(cleaned_line[1]), float(cleaned_line[2]), float(cleaned_line[3]))
    
    def clean_line(self, line):
        splitline = line.split(" ")
        # Spaces are removed, empty strings remain
        while("" in splitline): 
            splitline.remove("")
        return splitline
    
    def map_lines_to_elements(self,lines):
        elements = {}
        for line in lines:
                elements.update(self.map_line_to_element(line))
        return elements
    
    def map_line_to_element(self,line):
        splitline = line.split(",")
        return({splitline[2] : float(splitline[3])})
        
class Interpreter:
    
    def interpret(self, atoms, attachment_atom_name, crystal_atom = None, replacement_atom = None):
        if crystal_atom == None:
            self.atoms = atoms
            self.attachment_atom = self._find_molecule_attachment(self.atoms, attachment_atom_name)
            self._translate_molecule_to_origin(self.atoms, self.attachment_atom)
            self._remove_molecule_attachment(self.atoms, attachment_atom_name)
            self.center_of_mass_vector =  self._molecule_center_of_mass_vector(self.atoms, self.attachment_atom)

            output_content = [self.atoms, self.center_of_mass_vector]
        
        else:
            self.atoms = atoms
            self.attachment_atoms = self._find_crystal_ligand_positions(self.atoms, attachment_atom_name)
            self.center = self._find_crystal_center(self.atoms, crystal_atom)
            self.normal_vectors = self._find_crystal_normalized_normal_vectors(self.attachment_atoms, self.center)
            self._replace_crystal_attachment_atom(self.atoms, attachment_atom_name, replacement_atom)

            output_content = [self.atoms, self.attachment_atoms, self.normal_vectors]
        return output_content

    def _find_crystal_ligand_positions(self,atoms,attachment_atom):
        """ Find the atoms where ligands should be attached"""
        ligand_atoms = []
        for atom in atoms:
            if atom.kind == attachment_atom:
                ligand_atoms.append(atom)
        return ligand_atoms

    def _find_crystal_center(self,content,crystal_atom):
        """Finds the center of the crystal"""
        x_sum = 0
        y_sum = 0
        z_sum = 0
        crystal_atoms = 0

        for atom in content:
            if atom.kind == crystal_atom:
                crystal_atoms += 1
                x_sum += atom.x
                y_sum += atom.y
                z_sum += atom.z
        
        x_average = x_sum / crystal_atoms
        y_average = y_sum / crystal_atoms
        z_average = z_sum / crystal_atoms

        return Vector(x_average, y_average, z_average)

    def _find_crystal_normalized_normal_vectors(self, attachment_atoms, center):
        """Finds normal vectors for given positions and a center"""
        normal_vectors = []

        for attachment_atom in attachment_atoms:
            normal_vector = np.subtract([attachment_atom.x, attachment_atom.y, attachment_atom.z],center.vector)
            normal_vector = Vector(normal_vector[0], normal_vector[1], normal_vector[2])
            normal_vector.normalize()
            normal_vectors.append(normal_vector)
        return normal_vectors

    def _replace_crystal_attachment_atom(self, atoms, attachment_atom, replacement_atom):
        for atom in atoms:
            if atom.kind == attachment_atom:
                atom.kind = replacement_atom

    def _translate_molecule_to_origin(self, atoms, attachment_atom):
        """ Transforms ligand attachment site to origin """
        translate_to_origin_vector = Vector(- attachment_atom.x, - attachment_atom.y, - attachment_atom.z)
        mh._translate(atoms, translate_to_origin_vector)

    def _find_molecule_attachment(self,content,attachment_atom_name):
        for atom in content:
            if atom.kind == attachment_atom_name:
                return atom
        raise Exception("Ligand has no attachment site.")

    def _remove_molecule_attachment(self,atoms,attachment_atom_name):
        for atom in atoms:
            if atom.kind == attachment_atom_name:
                atoms.remove(atom)

    def _molecule_center_of_mass_vector(self,atoms,attachment_atom):
        read_csv = Reader(1)
        elements = read_csv.read(path.join('.','input','Periodic Table of Elements.csv'))
        center_of_mass_vector = Vector(0.,0.,0.)
        for atom in atoms:
            mass = elements[atom.kind]
            center_of_mass_vector.vector[0] += atom.x * mass
            center_of_mass_vector.vector[1] += atom.y * mass
            center_of_mass_vector.vector[2] += atom.z * mass

        #normalization
        center_of_mass_vector.normalize()
        return center_of_mass_vector

class Combiner:
    def combine(self, output):
        self.molecule_atoms = output[0][0]
        self.center_of_mass_vector = output[0][1]
        self.crystal_atoms = output[1][0]
        self.attachment_atoms = output[1][1]
        self.normal_vectors = output[1][2]
        self.molecules = self._make_molecules(self.molecule_atoms, self.center_of_mass_vector, self.attachment_atoms, self.normal_vectors)
        output_content = self._combine_molecules_with_crystal(self.molecules, self.crystal_atoms)
        return output_content

    def _make_molecules(self, molecule_atoms, center_of_mass_vector, attachment_atoms, normal_vectors):
        """Produces instances of molecules that can be attached to the crystal"""
        molecules = []
        for index, atom in enumerate(attachment_atoms):
            rotation_matrix = mh._make_rotation_matrix(center_of_mass_vector.vector, normal_vectors[index].vector)
            molecule = copy.deepcopy(molecule_atoms)
            rotated_molecule = mh._rotate(molecule, rotation_matrix)
            translated_molecule = mh._translate(rotated_molecule,Vector(atom.x,atom.y,atom.z))
            molecule = Molecule(translated_molecule)
            molecules.append(molecule)
        return molecules


    def _combine_molecules_with_crystal(self, molecules, crystal_atoms):
        crystal_with_attachments = []
        for molecule in self.molecules:
            for atom in molecule.atoms:
                crystal_with_attachments.append(atom)
        for atom in crystal_atoms:
            crystal_with_attachments.append(atom)
        return crystal_with_attachments
                   
class Writer:

    def __init__(self, delimiter):
        self.delimiter = delimiter

    def write(self,output, output_path):
        with open(output_path, 'w') as stream:
            stream.write(str(self.count_atoms(output))+ "\n\n")
            for atom in output:
                stream.write(atom.kind + self.delimiter + str(atom.x) + self.delimiter + str(atom.y) + self.delimiter + str(atom.z)+ "\n")
    
    def count_atoms(self,content):
        atoms = len(content)
        return atoms

def main(crystal_atom = "Si", attachment_atom_molecule = "P", attachment_atom_crystal = "H", replacement_atom = "O"):

    parser = Parser(Reader(2),Interpreter(), Combiner(),  Writer(' '))
    parser.interpret_molecule(glob.glob(path.join('input/molecule', '*.xyz'))[0],attachment_atom_molecule)
    parser.interpret_crystal(glob.glob(path.join('input/crystal', '*.xyz'))[0], attachment_atom_crystal,crystal_atom,replacement_atom)
    parser.combine()
    parser.write(path.join('output','output.xyz'))

if __name__ == '__main__':
    main()