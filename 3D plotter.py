# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 16:45:48 2023

@author: alexa
"""

from scipy.special import sph_harm
import numpy as np
from matplotlib import cm, colors
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


phi = np.linspace(0, np.pi, 250)
theta = np.linspace(0, 2*np.pi, 250)
phi, theta = np.meshgrid(phi, theta)



m, l = 0, 1

# Calculate the spherical harmonic Y(l,m) and normalize to [0,1]
fcolors = np.abs(sph_harm(m, l, theta, phi))
#fcolors_bar = np.conjugate(fcolors)
#fcolors = np.sqrt(np.real(fcolors*fcolors_bar))
fmax, fmin = fcolors.max(), fcolors.min()


x = fcolors*np.sin(phi) * np.cos(theta)
y = fcolors*np.sin(phi) * np.sin(theta)
z = fcolors*np.cos(phi)

fcolors = (fcolors - fmin)/(fmax - fmin)

fig = plt.figure()
ax = plt.axes(projection = '3d')
ax.view_init(20,30)
ax.plot_surface(x, y, z,  rstride=1, cstride=1, facecolors=cm.viridis(fcolors))
