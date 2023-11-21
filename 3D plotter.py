# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 16:45:48 2023

@author: alexa
"""

from scipy.special import sph_harm
import numpy as np
from matplotlib import cm, colors
import matplotlib.pyplot as plt
#from mpl_toolkits import mplot3d


phi = np.linspace(0, np.pi, 150)
theta = np.linspace(0, 2*np.pi, 150)
phi, theta = np.meshgrid(phi, theta)

Energies = [-13.612844276428223, -3.37939453125, -3.392089843749999]
ls = [0, 0, 1]

l = 1
ms = [-1,0,1]
fig = plt.figure()
pos =1

Coords = []
for m in ms:
    fcolors = np.abs(sph_harm(m, l, theta, phi))
    
    fmax, fmin = fcolors.max(), fcolors.min()
    
    
    X = fcolors*np.sin(phi) * np.cos(theta)
    Y = fcolors*np.sin(phi) * np.sin(theta)
    Z = fcolors*np.cos(phi)
    
    Coords = Coords + [(X,Y,Z)]
    
   # fcolors = (fcolors - fmin)/(fmax - fmin)

    # ax = fig.add_subplot( projection='3d')
    # pos = pos+1
    # ax.plot_surface(X, Y, Z, rstride=5, cstride=5, facecolors=cm.viridis(fcolors))
    
    # rotate the axes and update
    
for i in range(0,3):
    
for angle in range(0, 360):
    ax.view_init(20, angle)
    plt.draw()
    plt.pause(.01)