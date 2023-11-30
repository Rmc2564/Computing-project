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
from scipy.integrate import odeint, simpson

phi = np.linspace(0, np.pi, 150)
theta = np.linspace(0, 2*np.pi, 150)
phi, theta = np.meshgrid(phi, theta)

'The angular dependance is given by the spherical haronics, which are visualised'
'using this code'


l = 1
m = 1

sph = np.abs(sph_harm(m, l, theta, phi))
    
smax, smin = sph.max(), sph.min()
col = (sph - smin)/(smax - smin)

X = sph*np.sin(phi) * np.cos(theta)
Y = sph*np.sin(phi) * np.sin(theta)
Z = sph*np.cos(phi)

fig = plt.figure()

ax = fig.add_subplot( projection='3d')
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=cm.viridis(col))
ax.set_title('m = 1', fontsize = '15')