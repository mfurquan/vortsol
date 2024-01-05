#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 16:59:04 2023

@author: mfurquan
"""
# I3/(Re*I0) for a fixed cylinder

import numpy as np
from math import pi
import scipy.integrate as integ
import matplotlib.pyplot as plt

Re = 10.0
N  = 10000
delta = 1.0e-2
r = 0.5

def Vn(n):
    dth = 2.*pi/n
    eps = 2.0*r*np.sin(dth/2)
    def integrand1(ang):
        dr = np.sqrt(2.0*r*(r+delta)*(1.0-np.cos(ang))+delta**2)
        return np.cos(ang)*np.exp(-dr/eps)*r
    def integrand2(ang):
        n_dot_dr = (r+delta)*np.cos(ang) - r
        dr = np.sqrt(2.0*r*(r+delta)*(1.0-np.cos(ang))+delta**2)
        return n_dot_dr*(dr+eps)*np.exp(-dr/eps)*r/dr**2
    a = integ.quad(integrand1,0.,2.0*pi)
    b = integ.quad(integrand2,0.,2.0*pi)
    return a[0]/(Re*(2.0*pi*eps**2+eps*b[0]))

    # n_dot_I3 = 0
    # I0 = 2.0*pi*eps**2
    # for i in range(n):
    #     th = i*dth + dth/2
    #     dr = np.sqrt(2.0*r*(r+delta)*(1.0-np.cos(th))+delta**2)
    #     n_dot_dr = (r+delta)*np.cos(th) - r
    #     n_dot_I3 += np.cos(th)*np.exp(-dr/eps)*r*dth
    #     I0 += eps*n_dot_dr*(dr+eps)*np.exp(-dr/eps)*r*dth/dr**2
    #     return n_dot_I3/(Re*I0)
    
k = np.array(range(11,N+1))
v = np.empty([N-10])
for i in range(N-10):
    v[i] = Vn(i+11)
    
plt.plot(k,v)