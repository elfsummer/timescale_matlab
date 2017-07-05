"""
calcualte the adiabatic profile for Saturn
"""
"""
change log: 6/4/2015
add write_data()
"""

from matplotlib.pyplot import *
import sys
import numpy as np
from scipy import interpolate

def adiabat_saturn():

    # data for Cp
    T_data=np.array([100.,200.,250.,300.,350.,400.,450.,500.,])
    T_data = np.hstack([T_data,np.arange(600.,2200.,100.)])
    Cp_H2_data = np.array([28.154, 27.447, 28.344, 28.849, 29.081, 29.181, 29.229,
        29.260, 29.327, 29.441, 29.624, 29.881, 30.205, 30.581, 30.992, 31.423,
        31.861, 32.298, 32.725, 33.139, 33.537, 33.917, 34.280, 34.624]);
    Cp_He_data = 20.786
    X_H2 = 0.881
    Cp_data = Cp_H2_data*X_H2 + Cp_He_data*(1-X_H2)

    # temperature structure
    T_ref = np.hstack([134.8,np.arange(150.,2001.,25.)])
    T_mid = T_ref + 12.5

    # splin interpolation
    #Cp_fit=interpolate.splrep(T_data,Cp_data,s=0)
    #Cp=interpolate.splev(T_mid, Cp_fit, der=0)

    # cubic interpolation
    Cp_fit = interpolate.interp1d(T_data,Cp_data,kind='cubic')
    Cp = Cp_fit(T_mid)

    # pressure structure
    N = np.size(T_ref)
    log_T_ref = np.log10(T_ref)
    log_PG_ref = np.zeros(N)
    log_PG_ref[0] = 0.;

    R_gas = 8.31446
    for i in range(1,N):
        log_PG_ref[i]=log_PG_ref[i-1]+(log_T_ref[i]-log_T_ref[i-1])*Cp[i-1]/R_gas

    # specific heat interpolation
    # specific heat J/MOLE/K
    Cp_ref=Cp_fit(T_ref)
    # convert unit to J/KG/K
    mu=4*(1-X_H2)+2*X_H2
    Cp_ref=Cp_ref/(mu*1e-3)

    # thermal expansion coefficient
    alpha_ref = 1/T_ref

    # gravitational acceleration(cm/s^2)
    grav_acc_ref = 10.44*np.ones(N)

    # density, assuming ideal gas
    rho_ref = 10.**log_PG_ref*1e5*mu*1.e-3/(R_gas*T_ref)

    # pressure scale height
    H_p_ref = R_gas*T_ref/(mu*1.e-3*grav_acc_ref)

    # height informaiton
    z_ref = np.zeros(N)
    z_ref[0]=0.
    for i in range(1,N):
        z_ref[i]=z_ref[i-1] + Cp[i-1]/(mu*1e-3)/grav_acc_ref[i-1]*(
            T_ref[i]-T_ref[i-1])

    return T_ref,log_PG_ref,Cp_ref,alpha_ref,grav_acc_ref,rho_ref,H_p_ref,z_ref


"""generate the mesh for Saturn
"""
def mesh(N,T_0,T_N):

    # N is the numer of grids
    T_ref,log_PG_ref,Cp_ref,alpha_ref,grav_acc_ref,rho_ref,H_p_ref,z_ref = adiabat_saturn()
    
    # create mesh
    # find the index corrresponding to upper boundary 400 K
    index_0 = find_index(T_ref,T_0)
    index_N = find_index(T_ref,T_N)
    dz = (z_ref[index_N]-z_ref[index_0])/N
    z = dz/2 + dz*np.arange(N)+z_ref[index_0]

    # calculate the quantities at each grid
    T_fit = interpolate.interp1d(z_ref,T_ref,kind='cubic')
    T = T_fit(z)

    log_PG_fit = interpolate.interp1d(z_ref,log_PG_ref,kind='cubic')
    log_PG = log_PG_fit(z)
    # molecular weight for Saturn
    X_H2 = 0.881
    mu = 4.*(1-X_H2)+2.*X_H2
    # rho
    R_gas = 8.31446
    rho = 10**log_PG*1.e5*mu*1e-3/(R_gas*T)

    return z,dz,T,log_PG,rho,mu

def find_index(target_array,ref_element):
    """
    target_array is a sorted array with non-decreasing order
    find the index i in the target_array such that target_array[i-1] < ref_element <= target_array[i] 
    """
    i=0
    for element in target_array:
        if element >= ref_element:
            break
        i=i+1
    return i

#z,dz,T,log_PG,mu,rho = mesh(80)
#import matplotlib.pyplot as plt
#plt.clf()
#plt.plot(10**log_PG,T)
#plt.show()

def plot_mesh(N,T_0,T_N):
    z,dz,T,log_PG,rho,mu = mesh(N,T_0,T_N)
    figure()
    plot(10**log_PG,T)
    xlabel('PG')
    ylabel('T')
    show()
    return 0

def write_data():
    T_ref,log_PG_ref,Cp_ref,alpha_ref,grav_acc_ref,rho_ref,H_p_ref,z_ref = adiabat_saturn()
    PG_ref = 10**log_PG_ref
    import csv
    w = csv.writer(open("PT_Saturn.csv", "w"))
    for i in range(len(T_ref)):
        w.writerow([T_ref[i], PG_ref[i]])


