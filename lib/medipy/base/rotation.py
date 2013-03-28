##########################################################################
# MediPy - Copyright (C) Universite de Strasbourg
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

""" Conversions between different forms of 3D rotations

    * (axis, angle) where axis is normalized and angle measured in radians
    * rotation matrix, 3x3, orthogonal
    * normalized quaternion, (a, b, c, d) = a+b*i+c*j+d*k
"""

import math

import numpy

def axis_angle_to_matrix(axis, angle) :
    r""" Convert an (axis, angle) to a rotation matrix.
    
         This formula comes from Rodrigues' rotation formula,
         :math:`R = I + \hat{\omega} \sin \theta + \hat{\omega}^2 (1-\cos \theta)`
         where :math:`\hat{}` gives the antisymmetric matrix equivalent of the cross product
        
         .. math ::
            
             \hat{\omega} = \begin{matrix}
                                        0 & -\omega_z &  \omega_y \\
                                 \omega_z &         0 & -\omega_x \\
                                -\omega_y &  \omega_x &         0 \\
                            \end{matrix} 
        
         Diagonal terms can be rewritten :
         
         .. math ::
             
             \begin{matrix}
                 1+(1-\cos \theta)*(\omega_x^2-1) & = & 1+(1-\cos \theta)*\omega_x^2-(1-\cos \theta) \\
                                                  & = & \cos \theta+\omega_x^2*(1-\cos \theta)
             \end{matrix}
    """
    
    result = numpy.ndarray((3,3))
    
    cos = math.cos(angle)
    sin = math.sin(angle)
    one_minus_cos = 1.-cos
    
    result[0][0] = cos+axis[0]**2*(one_minus_cos)
    result[1][1] = cos+axis[1]**2*(one_minus_cos)
    result[2][2] = cos+axis[2]**2*(one_minus_cos)
    
    result[0][1] = -axis[2]*sin+axis[0]*axis[1]*one_minus_cos
    result[1][0] = +axis[2]*sin+axis[0]*axis[1]*one_minus_cos
    
    result[0][2] = +axis[1]*sin+axis[0]*axis[2]*one_minus_cos
    result[2][0] = -axis[1]*sin+axis[0]*axis[2]*one_minus_cos
    
    result[1][2] = -axis[0]*sin+axis[1]*axis[2]*one_minus_cos
    result[2][1] = +axis[0]*sin+axis[1]*axis[2]*one_minus_cos
    
    return result

def matrix_to_quaternion(matrix):
    """ Convert a rotation matrix to a unit quaternion.
    
        cf. http://arxiv.org/abs/math/0701759 , p.4
    """

    d0 = matrix[0,0]
    d1 = matrix[1,1]
    d2 = matrix[2,2]

    # Diagonal terms of the matrix yield 4*q_r^2, 4*q_i^2, 4*q_j^2, 4*q_k^2
    q_r_squared_4 = 1.0 + d0 + d1 + d2
    q_i_squared_4 = 1.0 + d0 - d1 - d2
    q_j_squared_4 = 1.0 - d0 + d1 - d2
    q_k_squared_4 = 1.0 - d0 - d1 + d2
    
    # Since we are going to divide by one of (q_r^2, q_i^2, q_j^2, q_k^2),
    # choose the largest one to avoid numerical errors
    max_value = max((q_r_squared_4, q_i_squared_4, 
                     q_j_squared_4, q_k_squared_4))
    
    if max_value == q_r_squared_4 :
        # Use q_r^2 to get other terms
        q_r = 0.5*math.sqrt(q_r_squared_4)
        factor = 2.*math.sqrt(q_r_squared_4) # i.e. 4*q_r
        
        q_i = (matrix[2,1]-matrix[1,2])/factor
        q_j = (matrix[0,2]-matrix[2,0])/factor
        q_k = (matrix[1,0]-matrix[0,1])/factor
    elif max_value == q_i_squared_4 :
        # Use q_i^2 to get other terms
        q_i = 0.5*math.sqrt(q_i_squared_4)
        factor = 2.*math.sqrt(q_i_squared_4) # i.e. 4*q_i
        
        q_r = (matrix[2,1]-matrix[1,2])/factor
        q_j = (matrix[1,0]+matrix[0,1])/factor
        q_k = (matrix[2,0]+matrix[0,2])/factor
    elif max_value == q_j_squared_4 :
        # Use q_j^2 to get other terms
        q_j = 0.5*math.sqrt(q_j_squared_4)
        factor = 2.*math.sqrt(q_j_squared_4) # i.e. 4*q_j
        
        q_r = (matrix[0,2]-matrix[2,0])/factor
        q_i = (matrix[1,0]+matrix[0,1])/factor
        q_k = (matrix[2,1]+matrix[1,2])/factor
    else : # max_value == q_j_squared_4
        # Use q_k^2 to get other terms
        q_k = 0.5*math.sqrt(q_k_squared_4)
        factor = 2.*math.sqrt(q_k_squared_4) # i.e. 4*q_k
        
        q_r = (matrix[1,0]-matrix[0,1])/factor
        q_i = (matrix[2,0]+matrix[0,2])/factor
        q_j = (matrix[2,1]+matrix[1,2])/factor
    
    return (q_r, q_i, q_j, q_k)

def axis_angle_to_quaternion(axis, angle):
    """ Convert an (axis, angle) to a unit quaternion.
    """
    
    result = numpy.asarray((math.cos(angle/2.), 0., 0., 0.))
    result[1:] = numpy.multiply(axis, math.sin(angle/2.))
    return result

def quaternion_to_axis_angle(quaternion):
    """ Convert a unit quaternion to an (axis, angle).
    """
    angle = 2.*math.acos(quaternion[0])
    axis = numpy.divide(quaternion[1:], math.sin(angle/2.))
    
    return axis, angle

def quaternion_to_matrix(quaternion) :
    r""" Combination of quaternion_to_axis_angle and axis_angle_to_matrix
         The following equalities are used to obtain this :
         
         theta = 2 acos(qr)
         <=> q_r = cos(theta/2)
         
         cos(theta) = 2*cos^2(theta/2)-1
                    = 2*q_r^2-1
         
         For the diagonal terms, we have :
         cos(theta)+omega_x^2*(1-cos(theta)) 
             = cos(theta)+q_i^2/sin^2(theta/2)*(1-cos(theta))
             = cos(theta)+q_i^2/sin^2(theta/2)*(2*sin^2(theta/2))
             = cos(theta)+2*q_i^2
             = 2*q_r^2-1+2*q_i^2 (q is unit-length, so 1=q_r^2+q_i^2+q_j^2+q_k^2)
             = q_r^2+q_i^2-(q_r^2+q_i^2+q_j^2+q_k^2-q_r-q_i)
             = q_r^2+q_i^2-q_j^2-q_k^2
         
         For the off_diagonal terms, we have :
         
         .. math ::
         
             \begin{matrix}
                 \omega_z*\sin \theta & = & q_k/\sin \frac{\theta}{2}*\sin \theta \\
                                      & = & q_k/\sin \frac{\theta}{2}*2*sin\frac{theta}{2}*\cos \frac{theta}{2} \\
                                      & = & 2*q_k*\cos \frac{theta}{2} \\
                                      & = & 2*q_k*q_r
             \end{matrix}
         and :
         
         .. math ::
         
             \begin{matrix}
                 \omega_x*\omega_y*(1-\cos \theta) & = & q_i*q_j/\sin^2 \frac{\theta}{2}*(1-\cos \theta) \\
                                                   & = & q_i*q_j/\sin^2 \frac{\theta}{2}*(2*\sin^2 \frac{\theta}{2}) \\
                                                   & = & 2*q_i*q_j
         Hence :
         
         .. math ::
         
             \omega_z*\sin \theta+\omega_x*\omega_y*(1-\cos theta) = 2*(q_i*q_j+q_r*q_k)
    """
    matrix = numpy.ndarray((3,3))
    matrix[0][0] = quaternion[0]**2+quaternion[1]**2-quaternion[2]**2-quaternion[3]**2
    matrix[1][1] = quaternion[0]**2-quaternion[1]**2+quaternion[2]**2-quaternion[3]**2
    matrix[2][2] = quaternion[0]**2-quaternion[1]**2-quaternion[2]**2+quaternion[3]**2
    matrix[0][1] = 2*(quaternion[1]*quaternion[2]-quaternion[0]*quaternion[3])
    matrix[0][2] = 2*(quaternion[3]*quaternion[1]+quaternion[0]*quaternion[2])
    matrix[1][2] = 2*(quaternion[2]*quaternion[3]-quaternion[0]*quaternion[1])
    matrix[1][0] = 2*(quaternion[1]*quaternion[2]+quaternion[0]*quaternion[3])
    matrix[2][0] = 2*(quaternion[3]*quaternion[1]-quaternion[0]*quaternion[2])
    matrix[2][1] = 2*(quaternion[2]*quaternion[3]+quaternion[0]*quaternion[1])
    
    return matrix
