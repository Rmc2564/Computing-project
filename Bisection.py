# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 15:20:59 2023

@author: alexa
"""



def bisection(f,r0,r1,tolerance): #applies the bisection method to a function f with initial guesses rlow and rhigh
    error = 'No sign change in given interval'
    if f(r0)*f(r1) > 0:
        return error

    if f(r0) > 0:
        ra = r0
        rb = r1
    else:
        ra = r1
        rb = r0
    while abs(ra - rb) > tolerance:
        m = (ra+rb)/2
        if f(m) > 0:
            ra = m
        else:
            rb = m
    return m

def test(x):
    return x**2 - 1

t = bisection(test,0.5,1.2,0.0001)
print(t)