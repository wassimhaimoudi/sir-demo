# -*- coding: utf-8 -*-
"""
Created on Wed May 11 20:59:51 2022

@author: Wassim
"""

import numpy as np
import matplotlib.pyplot as plt

#rungekutta4
def solution(f , ti , tf , n , y0):
    h = (tf - ti)/n 
    T = [ti]*n
    y = [0]*n
    y[0] = y0
    for i in range (n-1):
        k1 = f( y[i],T[i])
        k2 = f(y[i] + (h/2)*k1 , T[i] + (h/2))
        k3 = f(y[i] + (h/2)*k2 , T[i] + (h/2))    
        k4 = f(y[i] + h*k3 , T[i] + h)
        y[i + 1] = y[i] + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
        T[i + 1] = T[i] + h
    return y

#SIR
def SIR(Y,t):
    return np.array([-beta*Y[0]*Y[1], beta*Y[0]*Y[1] - gamma*Y[1], gamma*Y[1]])

#valeurs_initiales
I0 = 0.001             #population_infectés_initiale
S0 = 1. - I0           #population_susceptibles_initiale
R0 = 0                 #population_recouvris_initiale

njours = 100           #nombre_de_jours

beta = 1./6.           #taux_de_contacts_par_jours
gamma = 1./14.         #taux_de_retablissement_par_jours

ti , tf = 0 , 99       
Y0= np.array([S0 - I0 , I0 , R0])
T = np.linspace( ti , tf , njours  )
Sol = solution(SIR , ti , tf , njours, Y0) #solution_de_SIR
plt.plot( T , Sol)
plt.legend()
plt.grid(color = 'gray' , linestyle = '--' , linewidth = 0.5)
plt.title('Modèle SIR', fontdict=())
plt.xlabel('Nombres des jours')
plt.ylabel('Population')
plt.show()


    