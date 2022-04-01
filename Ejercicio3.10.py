# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 10:06:06 2022
@author: Eric-Portatil
Clase practica nanoelectronica 1/04/22
Sale lo que tendria que salir pero como si estuviera centrado en el 4/-4 preguntar en clase
"""
#%%
""""Importacion de librerias"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const 



#%%
"""Definición de constantes"""
kb=const.Boltzmann #J/K
eta=0.5 #adimensional
q=const.e #C
homo= -5.5 #eV
lumo=-1.5
hbar=const.hbar #J/S
T=298 #K
kb=kb/q #eV/K
hbar=hbar/q #eV/S


alfa=0.01 #adimensional

gammas=0.1 #eV
gammad=0.1 #eV
taus=hbar/gammas #S
taud=hbar/gammad #S
taut=taus+taud #S



#%%
"""Defino la función de Fermi"""
def Fermi(E,mu):
    return (1.0/(np.exp((E-mu)/(kb*T))+1.0))

"""Defino la funcion de mus y mud"""
def mus(V):
    return (Ef+eta*V)

def mud(V):
    return (Ef+(eta-1)*V)

#%%
"""Comienzo a calcular """

Vmin=-10. #Potencial incicial
Vmax=10. #Potencial final
Paso=0.05 #Paso entre cada potencial
Pasostot=int((Vmax-Vmin)/Paso)+1 #Pasos totales dados entre vmin y vmax
V=np.linspace(Vmin,Vmax,Pasostot) #genero un vector de pasostot numeros comprendidos entre vmax y vmin de pa

I1=np.zeros(Pasostot)
I2=np.zeros(Pasostot)
I3=np.zeros(Pasostot)


   
i=0
qcs=1.0
Ef= -2.5
while(i<Pasostot):

    N0H=2.0 #2 pq al principio esta ocupado por 2 electrones
    N0L=0.0
    
    U=0.0
    Uant=-1.0
    
    mu_s=mus(V[i])
    mu_d=mud(V[i])
    
    while((np.abs(Uant-U))>1e-6):
        NH =2.0*(taud*Fermi(U+homo,mu_s )+(taus*Fermi(U+homo, mu_d)))/taut
        NL =2.0*(taud*Fermi(U+lumo,mu_s )+(taus*Fermi(U+lumo, mu_d)))/taut
        Uant=U
        U=qcs*((NH-N0H)+(NL-N0L))
        U=Uant+alfa*(U-Uant)
    
    I1[i]=2.0*q*(Fermi(U+lumo, mu_s)-Fermi(U+lumo, mu_d)+(Fermi(U+homo, mu_s)-Fermi(U+homo, mu_d)))/taut
    I1[i]/=1e-6    
    
    i+=1
    
    
i=0
qcs=1.0
Ef= -3.5
while(i<Pasostot):

    N0H=2.0 #2 pq al principio esta ocupado por 2 electrones
    N0L=0.0
    
    U=0.0
    Uant=-1.0
    
    mu_s=mus(V[i])
    mu_d=mud(V[i])
    
    while((np.abs(Uant-U))>1e-6):
        NH =2.0*(taud*Fermi(U+homo,mu_s )+(taus*Fermi(U+homo, mu_d)))/taut
        NL =2.0*(taud*Fermi(U+lumo,mu_s )+(taus*Fermi(U+lumo, mu_d)))/taut
        Uant=U
        U=qcs*((NH-N0H)+(NL-N0L))
        U=Uant+alfa*(U-Uant)
    
    I2[i]=2.0*q*(Fermi(U+lumo, mu_s)-Fermi(U+lumo, mu_d)+(Fermi(U+homo, mu_s)-Fermi(U+homo, mu_d)))/taut
    I2[i]/=1e-6    
    
    i+=1
        
i=0
qcs=1.0
Ef= -5.0
while(i<Pasostot):

    N0H=2.0 #2 pq al principio esta ocupado por 2 electrones
    N0L=0.0
    
    U=0.0
    Uant=-1.0
    
    mu_s=mus(V[i])
    mu_d=mud(V[i])
    
    while((np.abs(Uant-U))>1e-6):
        NH =2.0*(taud*Fermi(U+homo,mu_s )+(taus*Fermi(U+homo, mu_d)))/taut
        NL =2.0*(taud*Fermi(U+lumo,mu_s )+(taus*Fermi(U+lumo, mu_d)))/taut
        Uant=U
        U=qcs*((NH-N0H)+(NL-N0L))
        U=Uant+alfa*(U-Uant)
    
    I3[i]=2.0*q*(Fermi(U+lumo, mu_s)-Fermi(U+lumo, mu_d)+(Fermi(U+homo, mu_s)-Fermi(U+homo, mu_d)))/taut
    I3[i]/=1e-6    
    
    i+=1



#%%
"""Represento mis vectores I frente a V"""
fig = plt.figure(1)        # Creo el objeto figura
plt.plot(V, I1, 'r',label="$E_F=-2.5\,eV$")             # Represento corriente frente a potencial
plt.plot(V, I2, 'b',label="$E_F=-3.5\,eV$")             # Represento corriente frente a potencial
plt.plot(V, I3, 'g',label="$E_F=-5.0\,eV$")             # Represento corriente frente a potencial
plt.xlabel('V (V)')       # Eje x
plt.ylabel('I (μA)')        # Eje y
plt.title('Corriente frente a diferencia de potencial')   # Título
plt.xlim([-10,10])
plt.legend()    

"""Preparo las variables para reprentar la conductancia"""
Vds=np.zeros(Pasostot-1)
G1=np.zeros(Pasostot-1)
G2=np.zeros(Pasostot-1)
G3=np.zeros(Pasostot-1)

i=1
while(i<(Pasostot-1)):
    Vds[i]=(V[i]-V[i-1])/2
    i1aux=(I1[i]-I1[i-1])
    i2aux=(I2[i]-I2[i-1])
    i3aux=(I3[i]-I3[i-1])
    
    G1[i]=i1aux/Vds[i]
    G2[i]=i2aux/Vds[i]
    G3[i]=i3aux/Vds[i]
    
    Vds[i]=V[i]
    i+=1

Vds[0]=Vmin

"""Represento la conductancia frente a V"""
fig = plt.figure(2)        # Creo el objeto figura
plt.plot(Vds, G1, 'r',label="$E_F=-2.5\,eV$")             # Represento corriente frente a potencial
plt.plot(Vds, G2, 'b',label="$E_F=-3.5\,eV$")             # Represento corriente frente a potencial
plt.plot(Vds, G3, 'g',label="$E_F=-5.0\,eV$")             # Represento corriente frente a potencial
plt.xlabel('V (V)')       # Eje x
plt.ylabel('G (μS)')        # Eje y
plt.title('Conductancia frente a diferencia de potencial')   # Título
plt.xlim([-10,10])
plt.legend()
