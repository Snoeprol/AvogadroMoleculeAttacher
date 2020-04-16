import os
import numpy as np
import sys
import unittest
from src.parser import *
from src import mathhelper


class TestGeometry(unittest.TestCase):
    # See if the right matrix is generated between two vectors
    def test_rotation(self):
        v1 = np.random.rand(3)
        v2 = np.random.rand(3)

        v1 = v1 / np.linalg.norm(v1)
        v2 = v2 / np.linalg.norm(v2)

        rotMatrix = mathhelper._make_rotation_matrix(v1,v2)
        v_1_rotated = np.dot(rotMatrix, v1)
        for i in range(len(v_1_rotated)):
            self.assertTrue(v_1_rotated[i] - v2[i] < 10E-5)
        self.assertTrue(np.linalg.det(rotMatrix) - 1 < 10E-5)

    def test_translation_translates_atoms_as_expected(self):
        # Arrange
        amount_of_atoms = np.random.randint(1,10000)
        generated_atoms = []
        expected_atoms_after_translate = []
        translation_vector = Vector(np.random.rand(), np.random.rand(), np.random.rand())
        for i in range(amount_of_atoms):
            scalar = np.random.randint(1,10)
            random_x = np.random.rand()* scalar
            random_y = np.random.rand()* scalar
            random_z = np.random.rand()* scalar

            random_atom = Atom(i, random_x, random_y, random_z)

            expected_x = random_x + translation_vector.get_x()
            expected_y = random_y + translation_vector.get_y()
            expected_z = random_z + translation_vector.get_z()

            expected_atom_after_translate = Atom(i, expected_x, expected_y, expected_z)

            generated_atoms.append(random_atom)
            expected_atoms_after_translate.append(expected_atom_after_translate)

        mathhelper._translate(generated_atoms, translation_vector)
        for index ,atom in enumerate(generated_atoms):

            assert (atom.kind == expected_atoms_after_translate[index].kind)
            assert (atom.x - expected_atoms_after_translate[index].x < 10E-8)
            assert (atom.y - expected_atoms_after_translate[index].y < 10E-8)
            assert (atom.z - expected_atoms_after_translate[index].z < 10E-8)
'''
class TestMolecularOperations(unittest.TestCase):
    def test_molecule_reader(self):

        reader = Reader('2')
        molecule = Interpreter(reader.read(os.path.join(r"D:\Downloads\GM\GM1\DFT\CryAttacher\input\molecule","pyrene.xyz")))

    def test_translateToOrigin_correctpath_writescontenttodisk(self):
        # Arrange
        input_path = "./tests/input/test.xyz"
        output_name = "./tests/output/pyreneOrigin.xyz"
        atom = 'P'
        with open('./tests/input/expected_output.xyz') as stream:
            expected_content = stream.readlines()

        # Act
        center.translateToOrigin(input_path, output_name, atom)
        
        # Assert
        with open(output_name, 'r') as stream:
            actual_content = stream.readlines()
        
        # Remove produced file 
        os.remove(output_name)

        assert(actual_content == expected_content)
    '''
if __name__ == '__main__':
    unittest.main()
#rotation_test()
#translation_translates_atoms_as_expected()
