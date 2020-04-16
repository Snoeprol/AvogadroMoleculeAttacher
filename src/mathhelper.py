import numpy as np

def _rotate(atoms, rotation_matrix):
    """ Rotate a set of atoms using a rotation_matrix""" 
    # implement rotate function from other file
    for atom in atoms:
        rotated_point = np.dot(rotation_matrix, [atom.x,atom.y,atom.z])
        atom.x = rotated_point[0]
        atom.y = rotated_point[1]
        atom.z = rotated_point[2]
    return atoms

def _make_rotation_matrix(vector_1,vector_2):
    """" Generates the rotation matrix from vector_1 to vector_2"""
    # Use formula for rotation matrix: R = I + A + A^2 * b
    # https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    v = np.cross(vector_1,vector_2)
    c = np.dot(vector_1,vector_2)
    s = np.linalg.norm(v)  
    b = (1-c)/s**2
    
    # Generate rotation matrix
    A = np.zeros((3,3))
    A[0][1] += -v[2]
    A[0][2] += v[1]
    A[1][2] += -v[0]
    A[1][0] += v[2]
    A[2][0] += -v[1]
    A[2][1] += v[0]
    B = np.dot(A,A)
    I = np.identity(3)

    R = I + A + b * B
    return(R)

def _translate(atoms, translation_vector):
    '''Translates content by vector'''
    for atom in atoms:
        atom.x += translation_vector.get_x()
        atom.y += translation_vector.get_y()
        atom.z += translation_vector.get_z()
    return atoms