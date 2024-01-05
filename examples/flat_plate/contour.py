#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 11:31:57 2023

@author: mfurquan
"""

import numpy as np
import matplotlib.pyplot as plt

N = 401
AR = 250.0
a = 0.5/AR
b = 0.5 - a

s = np.linspace(0, 4*b+2*np.pi*a,N)

def plate(s,a,b):
    if(s < np.pi*a/2):
        theta = s/a
        return np.array([b+a*np.cos(theta),a*np.sin(theta)])
    elif(s < 2*b+np.pi*a/2):
        return np.array([b+np.pi*a/2-s,a])
    elif(s < 2*b+3*np.pi*a/2):
        theta = (s-2*b)/a
        return np.array([-b+a*np.cos(theta),a*np.sin(theta)])
    elif(s < 4*b+3*np.pi*a/2):
        return np.array([s-3*b-3*np.pi*a/2,-a])
    else:
        theta = (s-4*b)/a
        return np.array([b+a*np.cos(theta),a*np.sin(theta)])

def trans(s,a,b):
    L = 4*b + 2*np.pi*a
    return (4*np.pi*s/L - np.sin(np.sin(np.sin(np.sin(np.sin(4*np.pi*s/L))))))*L/(4*np.pi)

r = np.empty([N,2])
for i in range(N):
    #r[i,:] = plate(s[i],a,b)
    r[i,:] = plate(trans(s[i],a,b),a,b)
    
y = np.empty([N])
for i in range(N):
    y[i] = trans(s[i],a,b)

plt.plot(s,y,'.-')
plt.show()
plt.plot(r[:,0],r[:,1],'.-')