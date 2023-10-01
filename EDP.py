#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TRABALHO FINAL 
MÉTODOS MATEMÁTICOS (PMC03)

Alunos:
    Gustavo Pereira
    Pedro Ignacio
    Rodrigo Kalko
    
Resolução da EDP
"""

import matplotlib.pyplot as plt
from scipy.integrate import simpson
import numpy as np

#Definição das constantes
a=1
b=1
step=0.01
q_ent=100
q_sai=-q_ent
k=10
NT=10

#Definindo os limites
x=np.arange(0, a+step, step)
y=np.arange(0, a+step, step)
x_axis, y_axis=np.meshgrid(x, y)

#Definindo Cn:
def Cn(n, a, b, q_ent, q_sai, k):
    x_primeira=np.arange(0, (a+step)/2, step)
    x_segunda=np.arange(a/2, a+step, step)
    y_primeira=[]
    y_segunda=[]
    for i in x_primeira:
        y_primeira.append(-q_ent/k*np.sin((n*np.pi*i)/a))
    for i in x_segunda:
        y_segunda.append(-q_sai/k*np.sin((n*np.pi*i)/a))
    Cn=2/(n*np.pi*np.cosh(n*np.pi*b/a))*(simpson(y_primeira, x_primeira)+simpson(y_segunda, x_segunda))
    return Cn

# Definindo a Temperatura:
T=np.zeros((len(x_axis), len(y_axis)))

for idi, i in enumerate(x):
    for idj, j in enumerate(y):
        for n in range(1, NT+1, 1):
            T[idj][idi]=T[idj][idi]+(Cn(n, a, b, q_ent, q_sai, k)*np.sin(n*np.pi*i/a)*np.sinh(n*np.pi*j/a))

fig,ax =plt.subplots(subplot_kw={"projection": "3d"})
surf=ax.plot_surface(x_axis, y_axis, T, cmap='plasma')
fig.colorbar(surf, pad=0.2)
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_zlabel('Temperature [°C]')
plt.show()