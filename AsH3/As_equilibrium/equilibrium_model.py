"""generic equilibrium model:
take a T-P profile and compute the equilibirum abundances along the profile
"""

import sys
import numpy as np

import cantera as ct

# read data from file
file_name = 'PT_Saturn.csv'
# SI units 
T, P = np.loadtxt(file_name,delimiter = ',', skiprows=0,usecols = (0,1),unpack = True)
# convert to SI units
P = P*1e5
print T,P

# input data
mechanism_directory = '/Users/wangdong/Documents/Cornell/2014_summer/cantera_input_files'
# mechanism name
#name = 'my_mechanism.cti'
name = 'As_thermo_data.cti'
# composition
# Jupiter AsH3 0.55 times enrichment
# composition = 'H2:1.0,He:0.157,CH4:2.37e-3,H2O:7.518e-3,NH3:6.64e-4,AsH3:2.2e-10'
# Saturn AsH3 7.5 times enrichment
composition = 'H2:1.0,He:0.135,CH4:5.33e-3,H2O:9.8e-3,NH3:4.54e-4,AsH3:3.0e-9'
# list of species that you want to plot
species_names = ['As','As2','As3','As4','AsH','AsH2','AsH3']
# output
output_file_name = 'As_eq_Sat.dat'

################################################################################
# no need to change the scripts below
# initialization the phase
ct.add_directory(mechanism_directory)
gas1=ct.Solution(name)
gas1.X = composition
nspecies = gas1.n_species
N = len(P)
x = np.zeros((nspecies,N))
# equilibrium states
for i in range(N):
    gas1.TP = T[i],P[i]
    gas1.equilibrate('TP')
    x[:,i]=gas1.X


# if plot argument exists, make a plot
if '--plot' in sys.argv[1:]:
    import matplotlib.pyplot as plt
    #species_names = ['CO','H2O','CH4','NH3','N2']
    #species_names = ['PH3','PH2','HOPO','H3PO4','PH','P2']
    n = len(species_names)
    species_indices = []
    for species_name in species_names:
        species_indices.append(gas1.species_index(species_name))

    plt.clf()
    for species_index, species_name in zip(species_indices,species_names):
        plt.plot(x[species_index,:],T,label = species_name)
    plt.xscale('log')
    #plt.yscale('log')
    plt.xlim(1e-20,1e-5)
    plt.gca().invert_yaxis()
    plt.xlabel('X')
    plt.ylabel('T(K)')
    plt.legend()
    plt.show()
    #plt.savefig('As_eq.png')
else:
    print("To view a plot of these results, run this script with the option --plot")
                
if '--output' in sys.argv[1:]:
    header_text = 'equilibrium mole fractions of Asenic containing species\n'
    header_variable_text = 'PG(pa), T(K), AsH3, AsH2'
    output_names = ['AsH3','AsH2']
    names_index = []
    output_data = np.zeros((len(output_names),N))
    for name in output_names:
        names_index.append(gas1.species_index(name))
    for i in range(len(output_names)):
        output_data[i,:]=x[names_index[i],:]
    output_data = np.concatenate((P[:,np.newaxis],T[:,np.newaxis],output_data.T),axis=1)
    np.savetxt(output_file_name,output_data,fmt='%.8e',header = header_text+header_variable_text)





    
        
    

