##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg, 2011-2012
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

"""
The rows of the matrices are the orientation of the other coordinate system in 
the LPS coordinate system. For example : 
  * the y axis of RAS is (0,-1,0), i.e. Anterior
  * the normal of the coronal slice (radiological convention) is (0,-1,0), 
    i.e. Anterior
  * the normal of the sagittal slice (radiological convention) is (0,0,-1), 
    i.e. Right

Since column-vector matrices would transform to LPS (cf. below) and since the
matrices are orthogonal (M^T=M^{-1}), the row-vector matrices transform from LPS
to the other coordinate system.

Let (e^1, e^2, e^3) and (e_1, e_2, e_3) be two bases of R^3. Let M be the matrix
formed by the column vectors of the coordinates of (e^1, e^2, e^3), expressed in 
(e_1, e_2, e_3), i.e.
    e^1_x e^2_x e^3_x 
M = e^1_y e^2_y e^3_y
    e^1_z e^2_z e^3_z

Let v be a column-vector :
    a
v = b
    c

By right-multiplying M by v, we get ae^1+be^2+ce^3, where the coordinates of the
e^i are expressed in (e_1, e_2, e_3). If the elements of v are coordinates 
expressed in (e^1, e^2, e^3), then M transforms a vector from (e^1, e^2, e^3)
to (e_1, e_2, e_3).

The preceding reasoning can be applied to row-vectors: let v_c and w_c be 
column-vectors, and let M be a matrix such that M.v_c = w_c. Let v_r (resp. w_r)
be the row vector corresponding to v_c (resp. w_c). Using the matrix-vector
multiplication formula, it is easy to prove that v_r.M = w_r.
"""

import numpy

# LPS coordinate system, in LPS coordinates. This is the DICOM coordinate system.
LPS = numpy.identity(3, dtype=numpy.float)
# RAS coordinate system, in LPS coordinates. This is the NIfTI coordinate system.
RAS = numpy.asarray([[1,  0,  0],
                     [0, -1,  0],
                     [0,  0, -1]], dtype=numpy.float)

# Transformation matrices from LPS to OpenGL coordinates, ZYX order
# The ASCII representations are what would be displayed on a screen.
slices = {
    "radiological" : {
        #      A
        #      ^
        #      |
        # R <--o--> L
        #      |
        #      v
        #      P
        "axial"    : [[-1,  0, 0],
                      [ 0, -1, 0],
                      [ 0,  0, 1]],
        #      S
        #      ^
        #      |
        # R <--o--> L
        #      |
        #      v
        #      I
        "coronal"  : [[0, -1, 0],
                      [1,  0, 0],
                      [0,  0, 1]],
        #      S
        #      ^
        #      |
        # A <--o--> P
        #      |
        #      v
        #      I
        "sagittal" : [[0, 0, -1],
                      [1, 0,  0],
                      [0, 1,  0]],
    },
    "neurological" : {
        #      A
        #      ^
        #      |
        # L <--o--> R
        #      |
        #      v
        #      P
        "axial"    : [[1,  0,  0],
                      [0, -1,  0],
                      [0,  0, -1]],
        #      S
        #      ^
        #      |
        # L <--o--> R
        #      |
        #      v
        #      I
        "coronal"  : [[0, 1,  0],
                      [1, 0,  0],
                      [0, 0, -1]],
        #      S
        #      ^
        #      |
        # P <--o--> A
        #      |
        #      v
        #      I
        "sagittal" : [[0,  0, 1],
                      [1,  0, 0],
                      [0, -1, 0]],
    },
}

def best_fitting_axes_aligned_matrix(direction):
    """ Compute the transformation matrix that best fits the given direction 
        matrix while being aligned to the axes.
    """
        
    transformation_matrix = numpy.zeros((3,3))
    for index, v in enumerate(direction) :
        max_alignment = None
        best_model = None
        for model in [(1,0,0), (0,1,0), (0,0,1)] :
            model = numpy.asarray(model)
            
            alignment = numpy.dot(v, model)
            
            if alignment > max_alignment :
                best_model = model
                max_alignment = alignment
            
            model = -model
            
            alignment = numpy.dot(v, model)
            
            if alignment > max_alignment :
                best_model = model
                max_alignment = alignment
        transformation_matrix[index,:] = best_model
    
    return transformation_matrix
