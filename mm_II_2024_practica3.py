import math

import numpy as np
from matplotlib import pyplot as plt

#Mecanica de Maquinaria
#Paralelo 104
#Grupo 1: Mecanismo b

def eclineal2(A,B):

    A = np.array(A)
    B = np.array(B)
    x = np.linalg.solve(A,B)
    return x



def nr2i(theta3k,theta4k,theta1,theta2,L1,L2,L3,L4,tol=0.01):
    error = tol+1
    while error>tol:
        F1 = L2*math.cos(theta2)+L3*math.cos(theta3k)-L1*math.cos(theta1)-L4*math.cos(theta4k)
        F2 = L2*math.sin(theta2)+L3*math.sin(theta3k)-L1*math.sin(theta1)-L4*math.sin(theta4k)
        J11 = -L3*math.sin(theta3k)
        J12 = +L4*math.sin(theta4k)
        J21 = +L3*math.cos(theta3k)
        J22 = -L4*math.cos(theta4k)
        d = eclineal2([[J11,J12],[J21,J22]],[-F1,-F2])
        theta3k += d[0]
        theta4k += d[1]
        error = math.sqrt((d[0]**2)+(d[1]**2))
    return theta3k,theta4k

N = 360
L1 = 6.00 #pulg #Longitud del eslabon 2
L2 = 3.00 #pulg
L3 = 10.00 #pulg
L4 = 12.00 #pulg

Rg2 = 1 #pulg
delta2 = 20 #grados
delta2 *= math.pi/180 #rad
Rg3 = 4 #pulg
delta3 = -30 #grados
delta3 *= math.pi/180 #rad
Rg4 = 6 #pulg
delta4 = 40 #grados
delta4 *= math.pi/180 #rad


theta1 = 0 #grados
theta1 *= math.pi/180 #rad
theta2 = np.linspace(2*math.pi,0,N+1) #rad
theta3 = np.zeros(N+1)
theta4 = np.zeros(N+1)

#Velocidad angular
omega2 = 20 #rpm
omega2 *= math.pi/30 #rad/s
tRev = 2*math.pi/abs(omega2) #tiempo de un revolucion [s]
t = np.linspace(0,tRev,N+1) #s
#Velocidad angular
omega3 = np.zeros(N+1)
omega4 = np.zeros(N+1)
#Aceleracion angular
#alpha2 es cero debido a que omega2 es constante
alpha3 = np.zeros(N+1)
alpha4 = np.zeros(N+1)

#Velocidad de las juntas
VA = np.zeros((N+1,3))
VB = np.zeros((N+1,3))
VC = np.zeros((N+1,3))
#Velocidad de los centros de masa
Vg2 = np.zeros((N+1,3))
Vg3 = np.zeros((N+1,3))
Vg4 = np.zeros((N+1,3))
#Aceleracion de las juntas
AA = np.zeros((N+1,3)) #Aceleracion junta del punto A
AB = np.zeros((N+1,3))
AC= np.zeros((N+1,3))
#Aceleracion de los centros de masa
Ag2 = np.zeros((N+1,3)) #Aceleracion del centro de masa del eslabon 2
Ag2mod = np.zeros((N+1))
Ag3 = np.zeros((N+1,3))
Ag3mod = np.zeros((N+1))
Ag4 = np.zeros((N+1,3))
Ag4mod = np.zeros((N+1))

for i in range(0,N+1):

    #Calculo de theta 3 y theta 4
    theta3[i], theta4[i] = nr2i(5.1,4.2,theta1,theta2[i],L1,L2,L3,L4)

    #Calculo de omega 3 y omega 4
    a11 = -L3*math.sin(theta3[i])
    a12 = +L4*math.sin(theta4[i])
    a21 = +L3*math.cos(theta3[i])
    a22 = -L4*math.cos(theta4[i])

    b1 = +L2*math.sin(theta2[i])*omega2
    b2 = -L2*math.cos(theta2[i])*omega2

    omega3[i],omega4[i] = eclineal2([[a11,a12],[a21,a22]],[b1,b2])

    #Calculo de alpha 3 y alpha 4
   
    b1 = L2*math.cos(theta2[i])*(omega2**2)+L3*math.cos(theta3[i])*(omega3[i]**2)-L4*math.cos(theta4[i])*(omega4[i]**2)
    b2 = L2*math.sin(theta2[i])*(omega2**2)+L3*math.sin(theta3[i])*(omega3[i]**2)-L4*math.sin(theta4[i])*(omega4[i]**2)
  
    alpha3[i],alpha4[i] = eclineal2([[a11,a12],[a21,a22]],[b1,b2])

    #Calculo de la velocidad lineal y centro de masa del eslabon 2
    omega2Vect = np.array([0,0,omega2])
    R2 = np.array([L2*math.cos(theta2[i]),L2*math.sin(theta2[i]),0])
    VA[i] = np.cross(omega2Vect,R2)
    R21 = np.array([Rg2*math.cos(theta2[i]+delta2),Rg2*math.sin(theta2[i]+delta2),0])
    Vg2[i] = np.cross(omega2Vect,R21)
    #Vcg4mod[i] = math.sqrt((Vcg2[i,0]**2)+(Vcg2[i,1]**2))

    #Calculo de la velocidad lineal y centro de masa del eslabon 3
    omega3Vect = np.array([0,0,omega3[i]])
    R3 = np.array([L3*math.cos(theta3[i]),L3*math.sin(theta3[i]),0])
    VB[i] = np.cross(omega3Vect,R3)
    R31 = np.array([Rg3*math.cos(theta3[i]+delta3),Rg3*math.sin(theta3[i]+delta3),0])
    Vg3[i] = np.cross(omega3Vect,R31)

    #Calculo de la velocidad lineal y centro de masa del eslabon 4
    omega4Vect = np.array([0,0,omega4[i]])
    R4 = np.array([L4*math.cos(theta4[i]),L3*math.sin(theta4[i]),0])
    VC[i] = np.cross(omega3Vect,R3)
    R41 = np.array([Rg4*math.cos(theta4[i]+delta4),Rg4*math.sin(theta4[i]+delta4),0])
    Vg4[i] = np.cross(omega4Vect,R41)

    #Calculo de la aceleracion lineal y centro de masa del eslabon 2
    AA[i] = np.cross(omega2Vect,VA[i])
    Ag2[i] = np.cross(omega2Vect,Vg2[i])
    Ag2mod[i] = math.sqrt((Ag2[i,0]**2)+(Ag2[i,1]**2))

    #Calculo de la aceleracion lineal y centro de masa del eslabon 3
    alpha3Vect = np.array([0,0,alpha3[i]])
    AB[i] = AA[i] + np.cross(alpha3Vect,R3) + np.cross(omega3Vect,VB[i])
    Ag3[i] = AA[i] + np.cross(alpha3Vect,R31) + np.cross(omega3Vect,Vg3[i])
    Ag3mod[i] = math.sqrt((Ag3[i,0]**2)+(Ag3[i,1]**2))

    #Calculo de la aceleracion lineal y centro de masa del eslabon 4
    alpha4Vect = np.array([0,0,alpha4[i]])
    AC[i] = AB[i] + np.cross(alpha4Vect,R4) + np.cross(omega4Vect,VC[i])
    Ag4[i] = AB[i] + np.cross(alpha4Vect,R41) + np.cross(omega4Vect,Vg4[i])
    Ag4mod[i] = math.sqrt((Ag4[i,0]**2)+(Ag4[i,1]**2))
  

plt.figure(1)
plt.subplot(1,2,1)
plt.plot(t,omega3)
plt.title('Velocidad angular del eslabon 3')
plt.xlabel('Tiempo [s]')
plt.ylabel('Omega 3 [rad/s]')
plt.grid(True)
plt.subplot(1,2,2)
plt.plot(t,alpha3)
plt.title('Aceleracion angular del eslabon 3')
plt.xlabel('Tiempo [s]')
plt.ylabel('Alpha 3 [rad/s^2]')
plt.grid(True)

plt.figure(2)
plt.subplot(1,2,1)
plt.plot(t,omega4)
plt.title('Velocidad angular del eslabon 4')
plt.xlabel('Tiempo [s]')
plt.ylabel('Omega 4 [rad/s]')
plt.grid(True)
plt.subplot(1,2,2)
plt.plot(t,alpha4)
plt.title('Aceleracion angular del eslabon 4')
plt.xlabel('Tiempo [s]')
plt.ylabel('Alpha 4 [rad/s^2]')
plt.grid(True)

plt.figure(3)
plt.subplot(1,3,1)
plt.plot(t,VA[:,[0,1]])
plt.title('Velocidad junta A')
plt.xlabel('Tiempo [s]')
plt.ylabel('Velocidad del punto A [pulg/s]')
plt.legend(['Velocidad en x','Velocidad en y'])
plt.grid(True)
plt.subplot(1,3,2)
plt.plot(t,VB[:,[0,1]])
plt.title('Velocidad junta B')
plt.xlabel('Tiempo [s]')
plt.ylabel('Velocidad del punto B [pulg/s]')
plt.legend(['Velocidad en x','Velocidad en y'])
plt.grid(True)
plt.subplot(1,3,3)
plt.plot(t,VC[:,[0,1]])
plt.title('Velocidad junta C')
plt.xlabel('Tiempo [s]')
plt.ylabel('Velocidad del punto C [pulg/s]')
plt.legend(['Velocidad en x','Velocidad en y'])
plt.grid(True)

plt.figure(4)
plt.subplot(1,3,1)
plt.plot(t,Vg2[:,[0,1]])
plt.title('Velocidad centro de masa eslabon 2')
plt.xlabel('Tiempo [s]')
plt.ylabel('Velocidad del punto g2 [pulg/s]')
plt.legend(['Velocidad en x','Velocidad en y'])
plt.grid(True)
plt.subplot(1,3,2)
plt.plot(t,Vg3[:,[0,1]])
plt.title('Velocidad centro de masa eslabon 3')
plt.xlabel('Tiempo [s]')
plt.ylabel('Velocidad del punto g3 [pulg/s]')
plt.legend(['Velocidad en x','Velocidad en y'])
plt.grid(True)
plt.subplot(1,3,3)
plt.plot(t,Vg4[:,[0,1]])
plt.title('Velocidad centro de masa eslabon 4')
plt.xlabel('Tiempo [s]')
plt.ylabel('Velocidad del punto g4 [pulg/s]')
plt.legend(['Velocidad en x','Velocidad en y'])
plt.grid(True)

plt.figure(5)
plt.subplot(1,3,1)
plt.plot(t,Ag2mod)
plt.title('Aceleracion del centro de masa eslabon 2')
plt.xlabel('tiempo [s]')
plt.ylabel('Magnitud aceleracion del punto g2 [in/s2]')
plt.grid(True)
plt.subplot(1,3,2)
plt.plot(t,Ag3mod)
plt.title('Aceleracion del centro de masa eslabon 3')
plt.xlabel('tiempo [s]')
plt.ylabel('Magnitud aceleracion del punto g3 [in/s2]')
plt.grid(True)
plt.subplot(1,3,3)
plt.plot(t,Ag4mod)
plt.title('Aceleracion del centro de masa eslabon 4')
plt.xlabel('Tiempo [s]')
plt.ylabel('Magnitud aceleracion del punto g4 [in/s2]')
plt.grid(True)

plt.figure(6)
plt.subplot(1,3,1)
plt.plot(t,AA[:,[0,1]])
plt.title('Aceleracion punto de union A')
plt.xlabel('Tiempo [s]')
plt.ylabel('Aceleracion punto A [pulg/s2]')
plt.legend(['Aceleracion en x','Aceleracion en y'])
plt.grid(True)
plt.subplot(1,3,2)
plt.plot(t,AB[:,[0,1]])
plt.title('Aceleracion punto de union B')
plt.xlabel('Tiempo [s]')
plt.ylabel('Aceleracion punto B [pulg/s2]')
plt.legend(['Aceleracion en x','Aceleracion en y'])
plt.grid(True)
plt.subplot(1,3,3)
plt.plot(t,AC[:,[0,1]])
plt.title('Aceleracion punto de union C')
plt.xlabel('Tiempo [s]')
plt.ylabel('Aceleracion punto C [pulg/s2]')
plt.legend(['Aceleracion en x','Aceleracion en y'])
plt.grid(True)
plt.show()