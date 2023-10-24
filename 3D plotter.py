# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 16:45:48 2023

@author: alexa
"""

from scipy.special import sph_harm
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


theta = np.linspace(0,np.pi,100000)
phi = np.linspace(0,2*np.pi,100000)
rs = np.linspace(0,2000,100000)

Y_20 = sph_harm(0,2, theta, phi)
Y_20bar = np.conj(Y_20)

def poltocart(r,theta,phi):
    x = r*np.cos(theta)*np.sin(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(phi)
    return [x, y, z]


new_coords = poltocart(rs, theta, phi)
x = new_coords[0]
y = new_coords[1]
z = new_coords[2]


fig = plt.figure()
ax = plt.axes(projection = '3d')
