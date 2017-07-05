# equilibrium plots for Ge contained species
# species interested:
# Ge, GeBr, GeBr2, GeCl, GeCl2, GeCl4, GeF, GeF2, GeH4, GeO, GeS, GeS2, Ge2, 
# Ge(cr), GeO2(II)

import matplotlib.pyplot as plt
import numpy as np

T = np.arange(300,1401,100)
N = T.size

T_GeF = T[3:]
GeF = np.array([0., 2.692e-19,3.807e-18, 1.809e-17, 6.090e-17, 1.581e-16, 3.350e-16, \
	6.046e-16, 9.596e-16])

T_GeF2 = T
GeF2 = np.array([2.425e-17, 2.555e-11, 4.888e-11, 1.602e-11, 7.215e-12, 9.155e-13, \
	1.143e-13,2.145e-14,5.348e-15,1.637e-15,5.832e-16,2.335e-16])

T_GeO2_II = T[0:2]
GeO2_II = np.array([7.6327e-9,7.6072e-9])

T_GeBr2 = T[0:9]
GeBr2 = np.array([0.,1.624e-20, 3.428e-19, 5.627e-19, 8.055e-19, 2.433e-19, 5.957e-20, \
	1.908e-20,0.])

T_GeCl = T[1:]
GeCl = np.array([0., 7.019e-20, 3.324e-17, 2.601e-15, 1.528e-14, 3.663e-14, 7.132e-14, \
1.181e-13, 1.720e-13, 2.258e-13, 2.725e-13])

T_GeCl2 = T[0:]
GeCl2 = np.array([0., 5.280e-17, 1.781e-13, 3.135e-13, 4.719e-13,  1.482e-13, 3.742e-14, 1.231e-14, \
4.847e-15, 2.167e-15, 1.061e-15, 5.571e-16])

T_GeCl4 = T[2]
GeCl4 = np.array([7.463e-19])

T_GeH4 = T[0:]
GeH4 = np.array([0., 3.535e-16, 1.417e-12, 7.913e-11, 1.4596e-9, 3.1161e-9, 3.6243e-9, 4.1870e-9, \
4.7628e-9, 5.3113e-9, 5.8008e-9, 6.2136e-9])

T_GeS = T[0:]
GeS = np.array([0.,7.333e-16, 3.288e-12, 1.714e-10, 2.6660e-9, 4.5152e-9, 4.0069e-9, 3.4415e-9,\
 2.8603e-9, 2.3035e-9, 1.8034e-9, 1.3786e-9])

T_GeO = T[1:]
GeO = np.array([0.,7.002e-20, 1.391e-16, 2.918e-14, 3.493e-13, 1.425e-12, 4.164e-12, 9.468e-12, \
1.770e-11, 2.836e-11, 4.016e-11])

T_Ge_cr = T[2:5]
Ge_cr = np.array([7.5790e-9,7.3659e-9,3.4994e-9])

T_GeS2 = T[3:]
GeS2 = np.array([0., 1.940e-19,1.035e-18,2.284e-18,4.120e-18,6.342e-18,8.594e-18,\
	1.050e-17,1.181e-17])

T_Ge = T[4:]
Ge = np.array([0., 1.101e-20,7.051e-19,1.866e-17, 2.573e-16,2.161e-15,1.235e-14, \
	5.202e-14])

T_GeBr = T[3:]
GeBr = np.array([0., 3.197e-19,2.447e-18,7.210e-18,1.655e-17,3.136e-17,5.109e-17,\
	7.372e-17,9.640e-17])

T_Ge2 = T[8:]
Ge2 = np.array([0., 4.346e-20,2.238e-19,8.310e-19])

# Ge, GeBr, GeBr2, GeCl, GeCl2, GeCl4, GeF, GeF2, GeH4, GeO, GeS, GeS2, Ge2, 
# Ge(cr), GeO2(II)

plt.close('all')
plt.figure()

plt.plot(Ge,T_Ge,label='Ge',lw = 2)
plt.xlim(1e-20,1e-8)
plt.ylim(300,1400)
plt.xscale('log')
plt.xlabel('mole fraction')
plt.ylabel('T')

plt.plot(GeBr,T_GeBr,label='GeBr',lw = 2)

plt.plot(GeBr2,T_GeBr2,label='GeBr2',lw = 2)

plt.plot(GeCl,T_GeCl,label='GeCl',lw = 2)

plt.plot(GeCl2,T_GeCl2,label='GeCl2',lw = 2)

#plt.plot(GeCl4,T_GeCl4,label='GeCl4')

plt.plot(GeF,T_GeF,label='GeF',lw = 2)

plt.plot(GeF2,T_GeF2,label='GeF2',lw = 2)

plt.plot(GeH4,T_GeH4,'--',label='GeH4',lw = 2)

plt.plot(GeO,T_GeO,'--',label='Geo',lw = 2)

plt.plot(GeS,T_GeS,'--',label='GeS',lw = 2)

plt.plot(GeS2,T_GeS2,'--',label='GeS2',lw = 2)

plt.plot(Ge2,T_Ge2,'--',label='Ge2',lw = 2)

plt.plot(Ge_cr,T_Ge_cr,'--',label='Ge(cr)',lw = 2)

plt.plot(GeO2_II,T_GeO2_II,'--',label='GeO2(II)',lw = 2)

plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

plt.show()


