# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 15:46:05 2023

@author: alexa
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def solve_coupled(sys, rs, y0, *args): #numerically solves a given differential equation
    sol = odeint(sys, y0, rs, *args)
    return sol


def test_sys(yv,t):
    x, y = yv
    dydt = [-5*x+4*y+np.exp(-2*t), -9*x+7*y+3*np.exp(-2*t)]
    return dydt

sol = solve_coupled(test_sys, np.linspace(0,5,1000), [0,0])

def y_anal(t):
    return 3*t*np.exp(t)

t = np.linspace(0,5,1000)
y_anal = y_anal(t)
y_err = ((sol[:,1] - y_anal)/y_anal)*100

plt.plot(np.linspace(0,5,1000), sol[:,0], label = 'x(t)')
plt.plot(np.linspace(0,5,1000), sol[:,1], label = 'y(t)')
plt.plot(t,y_anal)
plt.xlabel('time')
plt.legend()
plt.show()
plt.plot(t,y_err)
plt.show()