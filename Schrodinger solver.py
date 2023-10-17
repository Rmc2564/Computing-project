# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 15:20:59 2023

@author: alexa
"""

import numpy as np
from scipy.integrate import odeint, simpson
import matplotlib.pyplot as plt

alpha = 1/137
def coulomb(r):
    return -alpha/r

def schrodinger(y,r, E, l, mu):
    u,v = y
    dydr = [v, ((l*(l+1))/r**2)*u-2*mu*(E-coulomb(r))*u]
    return dydr

#Natural units in MeV

rs = np.linspace(1,2000,100000)
E = -1.360569*10**-5
l = 0
mu = 0.5110
y0 = [0,1]
a0 = 1/(mu*alpha)

solution = odeint(schrodinger, y0, rs, (E, l, mu), rtol = 1e-20) #Gives un-normalised wavefunction
wvfunc = solution[:,0]/rs
dist = rs**2*wvfunc*wvfunc
'now just need to normalise'

norm = simpson(dist, dx = 0.1)
pdist = dist/norm

plt.plot(rs, pdist)
E = E/4
l = 0
sol2 = odeint(schrodinger, y0, rs, (E,l,mu))
wvfnc2 = sol2[:,0]/rs
dist2 = rs**2*wvfnc2*wvfnc2

norm2 = simpson(dist2, dx = 0.1)
pdist2 = dist2/norm2

plt.plot(rs, pdist2)
