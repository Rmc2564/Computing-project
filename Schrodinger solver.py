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
    dydr = [v, (((l*(l+1))/r**2)-2*mu*(E-coulomb(r)))*u]
    return dydr

#Natural units in MeV

rs = np.linspace(3500,0.0000000000001,100000)
E = -1.304569*10**-5
l = 0
mu = 0.51099895000
y0 = [0,1]
a0 = 1/(mu*alpha)



solution = odeint(schrodinger, y0, rs, (E/4, l, mu), atol = 1e-20) #Gives un-normalised wavefunction
u = solution[:,0]
plt.plot(rs/a0,u)
plt.show()
norm = simpson(u*u, x = rs)
unorm = u/np.sqrt(abs(norm))

pdist = unorm*unorm
print(unorm[-1])
plt.plot(rs/a0,pdist)
plt.title('ground state')
plt.show()
'Now have a method for solving the schrodinger equation for a given energy,'
'just need to implement a modified bisection method'#


def getEnergy(E0, E1, tolerance, l):
    E0 = E0*10**-6
    E1 = E1*10**-6

    
    t1 = odeint(schrodinger, y0, rs, (E1, l, mu))[:,0]
    t0 = odeint(schrodinger, y0, rs, (E0, l, mu))[:,0]
    plt.plot(rs, t0, label = 'T0')
    plt.plot(rs, t1, label = 'T1')
    plt.legend()
    plt.show()
    norm1 = simpson(t1, x = rs)
    norm0 = simpson(t0, x = rs)
    t0 = t0/np.sqrt(norm0)
    t1 = t1/np.sqrt(norm1)
    test_0 = t0[-1]
    test_1 = t1[1]
    error = 'No sign change on interval'
    print('test_0: ' + str(test_0))
    print('test_1 : ' + str(test_1))
    if test_0*test_1 >0:
        return error
    
    if test_0 < 0:
        Elow = E0
        Ehigh = E1
    else:
        Ehigh  = E0
        Elow = E1
    Emid = (Ehigh+Elow)*0.5
    overflow_check = np.max(t1)
    print(abs(test_0 - test_1))
    while abs(test_0 - test_1) > tolerance:
        
        
        tlow = odeint(schrodinger, y0, rs, (Elow, l, mu))[:,0]
        normlow = simpson(tlow*tlow, x = rs[::-1])
        tlow = tlow/np.sqrt(normlow)
        test_0 = tlow[-1]

        thigh = odeint(schrodinger, y0, rs, (Ehigh, l, mu))[:,0]
        normhigh = simpson(thigh*thigh, x = rs[::-1])
        print(Emid*10**6)
        thigh = thigh/np.sqrt(normhigh)
        test_1 = thigh[-1]
        
       
        tmid = odeint(schrodinger, y0, rs, (Emid, l, mu))[:,0]
        normmid = simpson(tmid*tmid, x = rs[::-1])
        tmid = tmid/np.sqrt(normmid)
        
        testmid = tmid[-1]
        plt.plot(rs/a0, tmid)
        plt.title('Get smeared, idiot')
        if testmid < 0:
            Elow = Emid
        
        else:
            Ehigh = Emid
        Emid = (Ehigh + Elow)*0.5
    overflow_check = np.max(tmid)
    print(overflow_check)
    plt.show()
    if overflow_check > 1e7:
        return 'likely overflow'
    else:
        return Emid*10**6

    

E1 = getEnergy(-3.5,-3.8,0.0001, 0)
Sol_test = odeint(schrodinger, y0, rs, (E1*10**-6, l, mu))[:,0]



norm = simpson(Sol_test*Sol_test, x = rs[::-1])

dist2 = (Sol_test*Sol_test)/norm
plt.plot(rs/a0, dist2, label = 'numeric')

plt.title('n=2, l = 0')



