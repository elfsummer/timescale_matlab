% test 1: molecule GeH4
% results: pass
%{
% number of atoms
N = 5;
% temprature
T = 300;
% pressure
p = 1e5;
% parameters
mu = 76.64;
% parameters
% rotational constants
A = 79.9e9;
B = 79.9e9;
C = 79.9e9;
sigma = 12;
% parameters 
% vibrational frequencies
% length of mu = 3*N - 6 fo non-linear molecule, where N is the number of
% atoms, reference: Wang & Zhang, 2004
% unit: cm^-1
nu = [819 819 819 931 931 2106 2114 2114 2114];
% experimental results, consider scale factor = 1
sf_entropy = 1.0;
sf_cp = 1.0;

% number of unpaired electrons
n_unpaired = 0;
[S,Cp] = ...
    thermodynamic_computation(N,T,p,mu,A,B,C,sigma,nu,sf_entropy,sf_cp,...
    n_unpaired);

%}

% test 2 for radical: SiH2
% results: results are consistent with the nasa polynomial
% if we do not inlcude contribution from unpaired electrons, but why?

% number of atoms
N = 3;
% temprature
T = 300:100:1500;
% pressure
p = 1e5;
% parameters
mu = 30.1;
% parameters
% rotational constants (Hz)
A = 8.099*3e10;
B = 7.024*3e10;
C = 3.703*3e10;
sigma = 2;
% parameters 
% vibrational frequencies
% length of mu = 3*N - 6 fo non-linear molecule, where N is the number of
% atoms, reference: Wang & Zhang, 2004
% unit: cm^-1
nu = [999 1993 1996];
% experimental results, consider scale factor = 1
sf_entropy = 1.0;
sf_cp = 1.0;

% number of unpaired electrons
n_unpaired = 0;
for i=1:length(T)
[S(i),Cp(i)] = ...
    thermodynamic_computation(N,T(i),p,mu,A,B,C,sigma,nu,sf_entropy,sf_cp,...
    n_unpaired);
end

%{
% test on CH2: 
% tst results: also not including electronic partition function gives
% better results
% number of atoms
N = 3;
% temprature
T = 1000;
% pressure
p = 1e5;
% parameters
mu = 14.0;
% parameters
% rotational constants (Hz)
A = 19.805*3e10;
B = 11.249*3e10;
C = 7.239*3e10;
sigma = 2;
% parameters 
% vibrational frequencies
% length of mu = 3*N - 6 fo non-linear molecule, where N is the number of
% atoms, reference: Wang & Zhang, 2004
% unit: cm^-1
nu = [963 2806 3190];
% experimental results, consider scale factor = 1
sf_entropy = 1.0;
sf_cp = 1.0;

% number of unpaired electrons
n_unpaired = 0;
[S,Cp] = ...
    thermodynamic_computation(N,T,p,mu,A,B,C,sigma,nu,sf_entropy,sf_cp,...
    n_unpaired);
%}
