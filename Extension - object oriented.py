# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 21:07:57 2023

@author: alexa
"""
import numpy as np
from scipy.integrate import odeint, simpson
import matplotlib.pyplot as plt

rs = np.linspace(0.000000001, 4000, 10000)

alpha = 1/137
def coulomb(r):
        return -alpha/r

class radial:
    def __init__(self, initial, rs, energy, potential, l, mu):
        self.initial_conditions = initial
        self.energy = energy
        self.potential = potential
        self.rs = rs
        self.mass = mu
        self.l = l
    
    def schrodinger_solve(self):
        def schrodinger(y,r, E, l, mu):
            u,v = y
            dydr = [v, (((l*(l+1))/r**2)-2*mu*(E-self.potential(r)))*u]
            return dydr
    
        self.solution = odeint(schrodinger, self.initial_conditions, self.rs, (self.energy, self.l, self.mass))
        
        return self.solution
    

hydrogen_0 = radial([0,1], rs,-13.6*10**-6, coulomb, 0, 0.51099895000)

H0_solution = hydrogen_0.schrodinger_solve()[:,0]
H0_norm = np.abs(simpson((H0_solution), x = rs))
H0_solution = H0_solution/np.sqrt(H0_norm)

plt.plot(rs, H0_solution**2)