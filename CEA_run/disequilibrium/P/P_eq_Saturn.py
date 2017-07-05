"""Equilibrium Phosphorus model built for Saturn"""
import sys
import numpy as np

import cantera as ct

from mesh_saturn import *

from matplotlib.pyplot import *

# create mesh for Saturn
N = 70
[z,dz,T,log_PG,mu_average,rho] = mesh(N)
PG = 10**log_PG

# initialization the phase
ct.add_directory('/Users/wangdong/Documents/2014_summer/cantera_input_files')
gas1=ct.Solution('P_Reaction_set.cti')
nspecies = gas1.n_species
    
# O/H2 protosolar: 1.074e-3, 10 times protosolar 
#gas1.X = 'H2:1.0,He:0.135,H2O:1.074e-2,PH3:7.28e-6'
# O/H2 protosolar:  1.074e-3 20 times protosolar
gas1.X = 'H2:1.0,He:0.135,H2O:2.148e-2,PH3:7.28e-6'

y_eq = np.zeros((nspecies,N))
# set up the initial state using equilibrium states
for i in range(N):
    gas1.TP = T[i],PG[i]*1.e5
    gas1.equilibrate('TP')
    y_eq[:,i]=gas1.Y

# change mass fractions to mole fractions
x_eq = np.zeros((nspecies,N))
for i in range(N):
    x_eq[:,i] = y_eq[:,i]/gas1.molecular_weights*gas1.mean_molecular_weight

# make some plots
#option 1 
plot(x_eq[0,:],T,)
plot(x_eq[1,:],T,'g')
plot(x_eq[2,:],T,'r')
plot(x_eq[3,:],T,'k')
plot(x_eq[4,:],T,'y')
plot(x_eq[5,:],T,'b-')
#option 2
#for i in range(nspecies):
#    plot(x_eq[i,:],T,'--')

xscale('log')
xlabel('mole fraction')
ylabel('T')
xlim([1e-14,1e-5])
#legend('PH3','H3PO4(l)','H3PO4(cr)')
show()
#savefig('P_species_Jupiter')

# save arrays into file
np.savetxt('X_eq_with_solids.txt',x_eq)
np.savetxt('temp_20.txt',T)



