import numpy as np
import matplotlib.pyplot as plt

import openmdao.api as om

from openaerostruct.geometry.utils import generate_mesh
from openaerostruct.aerodynamics.eval_mtx import _compute_finite_vortex
from openaerostruct._SAFA_TEST.VORTEX import VORTEX
from openaerostruct._SAFA_TEST.VORTEX import VORTEX2

# Create a dictionary to store options about the mesh
mesh_dict = {"num_y": 3, "num_x": 2, "wing_type": "rect", "symmetry": True,"span": 1.0,"root_chord": 1.,}

# Generate the aerodynamic mesh based on the previous dictionary
mesh = generate_mesh(mesh_dict)

#print(mesh)


def plot_mesh(mesh):
    """ this function plots to plot the mesh """
    mesh_x = mesh[:, :, 0]
    mesh_y = mesh[:, :, 1]
    plt.figure(figsize=(6, 3))
    color = 'k'
    #plt.plot(np.linspace(-30,30,61),np.zeros(61))
    for i in range(mesh_x.shape[0]):
        plt.plot(mesh_y[i, :], mesh_x[i, :], color, lw=1)
        plt.plot(-mesh_y[i, :], mesh_x[i, :], color, lw=1)   # plots the other side of symmetric wing
    for j in range(mesh_x.shape[1]):
        plt.plot(mesh_y[:, j], mesh_x[:, j], color, lw=1)
        plt.plot(-mesh_y[:, j], mesh_x[:, j], color, lw=1)   # plots the other side of symmetric wing
    plt.axis('equal')
    plt.xlabel('span, m')
    plt.ylabel('chord, m')
    plt.legend()
    plt.show()



#plot_mesh(mesh)

'''
#OAS Core Model Test


#Distance Test
dists = np.linspace(1e-1,1e-12,101)
resultsOAS = []
resultsKfid = []
test = 1e-10


for test in dists:
    r1 = np.array([0.5,0,test])
    r2 = np.array([-0.5,0,test]) 

    result = _compute_finite_vortex(r1,r2)
    #print('OAS MODEL',result)
    resultsOAS.append(result[1])

    result = VORTEX2(r1,r2)
    #print('Kfid + Safa Model',result)
    resultsKfid.append(result[1])


plt.figure()
plt.plot(dists,resultsOAS,label='OAS')
plt.plot(dists,resultsKfid,label='OAS + Improved Core Model')
plt.xlabel('r(meters)')
plt.ylabel('Induced velocity(m/s)')
plt.xscale('log')
plt.grid()
plt.legend()
plt.savefig('CoreModelComparisonTol4.pdf')
plt.show()

'''

