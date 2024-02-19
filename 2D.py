# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 01:27:15 2019

@author: Pichau
"""

import numpy as np 
import seaborn 
seaborn.set(style='ticks')
import matplotlib.pyplot as plt
from random import *

M = 2
N= 120
R= 0.75
m = M/N
h= 0.04/np.sqrt(N/1000)
dt = 0.04
T = 400
ni = 1
k = 0.1
g = (8*k*M)/(np.power(R,4)*np.pi)



def r_0():
    r = []
    for i in range(N):
        r.append([uniform(-0.75,0.75),uniform(-0.75,0.75)])
    r = np.array(r)
    return r

    
def V_0():
    V = []
    for i in range(N):
        V.append([0,0])
    V = np.array(V)
    return V    

    

def ker(r,h):
    C = 5/(14*np.pi*np.power(h,2))
    q = np.linalg.norm(r)/h
    if 0<= q < 1:
        ker = C*(np.power(2-q,3) - 4*np.power(1-q,3))
        return ker
    if 1<= q < 2:
        ker = C*np.power(2-q,3)
        return ker
    if 2<= q:
        ker = 0
        return ker
    

def gradker(r,h):
    C = 5/(14*np.pi*np.power(h,2))
    norm = np.linalg.norm(r)
    q = norm/h
    if 0<= q < 1:
        ker = np.array([((3*C)/np.power(h,3)) * r[0]*(3*norm-4*h),((3*C)/(np.power(h,3))*r[1]*(3*norm-4*h))])
        return ker
    if 1<= q < 2:
        ker = np.array([((-3*C)/(np.power(h,3)*norm)) * r[0]*np.power((norm-2*h),2),((-3*C)/(np.power(h,3)*norm)) * r[1]*np.power((norm-2*h),2)])
        return ker
    if 2<= q:
        ker = 0
        return ker



def calc_rho(r):
    RHO =[]
    for i in range(N):
        RHO.append(m*ker(0,h))
    for i in range(N):
        
        for j in range(i+1,N):
            dij = r[i]-r[j]
            rhoij = m*ker(dij,h)
            RHO[i] += rhoij
            RHO[j] += rhoij
    return np.array(RHO)
  

def calc_a(r,V,RHO,P):
    A=[]
    for i in range(N):
        A.append([0,0])
    
    for i in range(N):
        A[i] += -ni*V[i] -g*r[i]
    for i in range(N):
        for j in range(i+1,N):
            dij = r[i]-r[j]
            aij = -m*((P[i]/(np.power(RHO[i],2)))+(P[j]/(np.power(RHO[j],2))))*gradker(dij,h)
            A[i] += aij
            A[j] += -aij
    return np.array(A)


def Graf(r,RHO):
    q = np.linspace(0.0, 0.75, 100)
    U = []
    rr = []
    X=[]
    Y=[]
    for i in range(len(q)):
        u = g/(4*k)*(np.power(R,2)-np.power(q[i],2))
        U.append(u)
    for i in range(N):
        a = np.linalg.norm(r[i])
        rr.append(a)
        Y.append(r[i][1])
        X.append(r[i][0])
    plt.figure(2)
    plt.plot(q,U,color = "cornflowerblue",label = "Teórico")
    plt.scatter(rr,RHO, color = "darkorange", label = "Simulado")
    plt.xlabel("Posição")
    plt.ylabel("Densidade")
    plt.title("t = {}".format(T))

    plt.grid(True, which='both')
    plt.legend()
    seaborn.despine(offset=0)
    
    
    plt.figure(3)
    plt.scatter(X,Y, c=RHO,cmap ="jet")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("t = {}".format(T))
    cbar = plt.colorbar()
    cbar.set_label('Densidade')
    plt.grid(True, which='both')
    plt.legend()
    seaborn.despine(offset=0)
    
    
    plt.show()

def main():
    r = r_0()
    V = V_0()
    RHO = calc_rho(r)
    
    P = k*np.power(RHO,2)
    A = calc_a(r,V,RHO,P)
    V_mhalf = V - 0.5*dt*A
    for t in range(T):
        V_phalf = V_mhalf +A*dt
        r += V_phalf*dt
        V = 0.5*(V_phalf+V_mhalf)
        V_mhalf = V_phalf
        RHO = calc_rho(r)
        P = k*np.power(RHO,2)
        A = calc_a(r,V,RHO,P)
        
        print(t)
    Graf(r,RHO)


main()