import numpy as np
import matplotlib.pyplot as plt

import openmdao.api as om

from openaerostruct.geometry.utils import generate_mesh


# Create a dictionary to store options about the mesh
mesh_dict = {"num_y": 3, "num_x": 2, "wing_type": "rect", "symmetry": True,"span": 1.0,"root_chord": 1.,}

# Generate the aerodynamic mesh based on the previous dictionary
mesh = generate_mesh(mesh_dict)

print(mesh)


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

plot_mesh(mesh)