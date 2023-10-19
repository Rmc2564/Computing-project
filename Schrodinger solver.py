# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 15:20:59 2023

@author: alexa
"""

import numpy as np
from scipy.integrate import odeint, simpson, solve_ivp
import matplotlib.pyplot as plt

alpha = 1/137
def coulomb(r):
        return -alpha/r

def schrodinger(y,r, E, l, mu):
    u,v = y
    dydr = [v, (((l*(l+1))/r**2)-2*mu*(E-coulomb(r)))*u]
    return dydr

#Natural units in MeV

rs = np.linspace(4000,1,100000)
E = -1.360569*10**-5
l = 0
mu = 0.5110
y0 = [0,1]
a0 = 1/(mu*alpha)



solution = odeint(schrodinger, y0, rs, (E, l, mu), atol = 1e-20) #Gives un-normalised wavefunction
u = solution[:,0][::-1]
rs = rs[::-1]
norm = simpson(u*u, x = rs)
unorm = u/np.sqrt(abs(norm))

pdist = unorm*unorm

plt.plot(rs/a0,pdist)

E = E/4
l = 0

solution_2_0 = odeint(schrodinger, y0, rs, (E, l, mu))
u20 = solution_2_0[:,0][::-1]
rs = rs[::-1]
norm20 = simpson(u20*u20, x = rs)
u20norm = u20/np.sqrt(abs(norm20))

pdist20 = u20norm*u20norm

plt.plot(rs/a0, pdist20)

