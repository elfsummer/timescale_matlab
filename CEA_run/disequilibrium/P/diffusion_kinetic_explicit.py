"""Diffsuion kinetic model built for Saturn"""
import sys
import numpy as np

import cantera as ct

from mesh_saturn import *

# create mesh for Saturn
N = 50
z,dz,T,log_PG,rho,mu = mesh(N,400,1100)
PG = 10**log_PG

# initialization the phase
ct.add_directory('/Users/wangdong/Documents/2014_summer/cantera_input_files')
gas1=ct.Solution('P_Reaction_set.cti')
nspecies = gas1.n_species

# restart or not?
restart = 0

if restart == 0:
    # set up initial run
    # O/H2 protosolar: 1.074e-3, 10 times protosolar 
    #gas1.X = 'H2:1.0,He:0.135,H2O:1.074e-2,PH3:7.28e-6'
    # O/H2 protosolar:  1.074e-3 10 times protosolar
    gas1.X = 'H2:1.0,He:0.135,H2O:1.074e-2,PH3:7.28e-6'
    y = np.zeros((nspecies,N))
    # set up the initial state using equilibrium states
    for i in range(N):
        gas1.TP = T[i],PG[i]*1.e5
        gas1.equilibrate('TP')
        y[:,i]=gas1.Y
    # record equilibrium abundances
    y_eq = y.copy()
    t_0 = 0
    # write equilibrium values into files
    np.savetxt('y_equilibrium.txt',y_eq)
else:
    # restart
    # read the starting time from file
    f=open('time.txt','r')
    t_0=int(f.read())
    f.close()
    # read the array from file
    y = np.loadtxt('y_t.txt')
    y_eq = np.loadtxt('y_equilibrium.txt')

print 'end of initilization'

# molecular weight (conver to kg/mole)
mu = gas1.molecular_weights*1.e-3

# eddy diffsuion coefficient
K = 1.e4*np.ones(N)

# diffusion initialization
dt = 1000
Tstep = 1000
# calculate some useful inermediates
rho_K_b = (rho[0:N-1]*K[0:N-1] + rho[1:N]*K[1:N])/2
y_new = np.zeros((nspecies,N))

# chemcial evolution initilzation
r1 = ct.IdealGasReactor(contents = gas1, energy = 'off')
network1=ct.ReactorNet([r1])
#network1.atol = 1e-10
#network1.rtol = 1e-3
network1.max_err_test_fails = 10

# time evolution starts here
for j in range(Tstep):
    # diffusion for dt/2
    for k in range(nspecies):
        # fisrt half time step: explicit
        for i in range(N):
            if i==0:
                phi_l = 0
                phi_u = rho_K_b[i]*(y[k,i+1]-y[k,i])/dz
            elif i==N-1:
                phi_u = 0
                phi_l = rho_K_b[i-1]*(y[k,i]-y[k,i-1])/dz;
            else:
                phi_u = rho_K_b[i]*(y[k,i+1]-y[k,i])/dz
                phi_l = rho_K_b[i-1]*(y[k,i]-y[k,i-1])/dz;
            y_new[k,i]= y[k,i] + (phi_u-phi_l)*dt/(2*dz*rho[i])
    
    y = y_new*1.0
    
    # chemical evolution for dt
    # number of sub iteration
   
    for i in range(N):
        gas1.TPY = T[i],PG[i]*1.e5,y[:,i]
        r1.syncState()
        network1.reinitialize()
        network1.set_initial_time(0.0)
        network1.advance(dt)
        y[:,i]=gas1.Y

    # diffusion for another dt/2
    for k in range(nspecies):
        # fisrt half time step: explicit
        for i in range(N):
            if i==0:
                phi_l = 0
                phi_u = rho_K_b[i]*(y[k,i+1]-y[k,i])/dz
            elif i==N-1:
                phi_u = 0
                phi_l = rho_K_b[i-1]*(y[k,i]-y[k,i-1])/dz;
            else:
                phi_u = rho_K_b[i]*(y[k,i+1]-y[k,i])/dz
                phi_l = rho_K_b[i-1]*(y[k,i]-y[k,i-1])/dz;
            y_new[k,i]= y[k,i] + (phi_u-phi_l)*dt/(2*dz*rho[i])
    
    y = y_new*1.0
    
    print j

    # record data for every 100 time steps
    if (j+1)%100 == 0:
        # save data for plot or restart
        t = t_0 + (j+1)*dt
        # save the time to file
        f = open('time.txt', 'w')
        f.write(str(t))
        f.close()
        # save arrays into file
        np.savetxt('y_t.txt',y)
    

# if plot argument exists, make a plot
if '--plot' in sys.argv[1:]:
    import matplotlib.pyplot as plt
    
    plt.clf()
    plt.plot(y[0,:],T,y[1,:],T,y[2,:],T)
    plt.xscale('log')
    plt.xlabel('mass fraction')
    plt.ylabel('T')
    plt.show()
else:
    print("To view a plot of these results, run this script with the option --plot")
                



    
        
    

