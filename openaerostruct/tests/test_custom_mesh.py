# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 14:27:47 2021

@author: rcharayr
"""

import numpy as np
import matplotlib.pyplot as plt
from openaerostruct.geometry.utils import gen_custom_mesh

"""A function to plot a mesh"""
def plotmesh(mesh):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_wireframe(mesh[:,:,0], mesh[:,:,1], mesh[:,:,2],rstride=1,cstride=1, label = 'initial') 
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.legend(loc=0)
    plt.show()

# rmq : soit on donne le taper ratio en entrée et on la chord_distrib est de dim(num_y-2) soit on ne donne pas le taper ratio en entrée et on utilise une chord_distrib de dim(num_y).  
# mesh = gen_custom_mesh(num_x = 5, num_y = 7, span = 60, chord_distrib = [6, 9, 13, 20, 13, 9, 6], sweep_angle = np.pi/6, dihedral_angle = 0.1, taper_ratio = 0.3, wing_twist_distrib = [0.1, 0.05, 0.01, 0, 0.01, 0.05, 0.1], span_cos_spacing= 0.5, chord_cos_spacing= 0.5)

# print(mesh)
# plotmesh(mesh)

from openaerostruct.geometry.utils import generate_mesh
#création d'un dictionnaires qui contient les options du maillage :
mesh_dict = {'num_y' : 7,
             'num_x' : 5,
             'wing_type' : 'customized',
             'symmetry' : False,
             'chord_distrib' : [6, 9, 13, 20, 13, 9, 6],
             'sweep' : np.pi/6,
             'dihedral_angle_distrib' : [0.1, 0.04, 0.02, 0.02, 0.04, 0.1],
             'taper' : 0.3,
             'wing_twist_distrib' : 4*[0.1, 0.05, 0.01, 0, 0.01, 0.05, 0.1],
             'chord_cos_spacing' : 1,
             'span_cos_spacing' : -1}

# génération d'un maillage basé sur le dictionnaire grâce à la fonction generate_mesh() :
mesh = generate_mesh(mesh_dict)

print(mesh)
plotmesh(mesh)