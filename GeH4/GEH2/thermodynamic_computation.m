function [S_total,Cp_total] = ...
    thermodynamic_computation(N,T,p,mu,A,B,C,sigma,nu,sf_entropy,sf_cp,n_unpaired)
% created 08/25/2015
% This function is to evalute thermodynamic parameters S and Cp
% from molecular parameters for non-linear molecules only. scripts for 
% linear molecules will be implemented later

%% input parameters
% number of atoms in the molecule
% T: temprature for the evaluation (K)
% p: pressure for the evaluation (Pa)
% mu: molecular weight (no unit or g/mol)
% A, B, C: rotational constants (Hz)
% sigma: rotational symmetry number
% nu: vibrational frequencies
% sf_entropy: scale factor for vibrational frequencies for entropy 
% sf_cp: scale factor for vibrational frequencies for molar heat capacity

%% output parameters
% molar entropy: S_total (J/mol/K)
% molar heat capacity: Cp_total (J/mol/K)


%% check if the number of frequencies are correct
if length(nu)~=3*N-6
    error('the number of vibrational frequencies are not correct!');
end

%% constants used
k = 1.38066e-23;
Na = 6.022137e23;
R = k*Na;
h = 6.626076e-34;
% speed of light in the unit of cm/s
c = 29979245800;
m_p = 1.67e-27;

%% translational partition function
m = mu*m_p;
% computation
S_trans = R*(3./2*log(2*pi*m/h^2) + 5./2*log(k*T) -log(p) +5./2);
Cp_trans = 5./2*R;

%% rotational partition function
% computation for non-linear molecules
S_rot = R*(3/2*log(k*T/h) - 0.5*log(A*B*C/pi)-log(sigma)+3./2);
Cp_rot = 3/2*R;

%% vibrational partion function
% convert to units of Hz
nu = nu*c;
% computation
nu_entropy = nu*sf_entropy;
nu_cp = nu*sf_cp;
n_freq = length(nu);
S_vib = 0;
Cp_vib = 0;
for i=1:n_freq
    S_vib = S_vib - R*log(1-exp(-h*nu_entropy(i)/(k*T))) + ...
        R*h*nu_entropy(i)/(k*T)*exp(-h*nu_entropy(i)/(k*T))...
        /(1-exp(-h*nu_entropy(i)/(k*T)));
    Cp_vib = Cp_vib + R*(h*nu_cp(i)/(k*T))^2*exp(-h*nu_cp(i)/(k*T))/...
        (1-exp(-h*nu_cp(i)/(k*T)))^2;
end


%% do not include electronic partation function so far! 
% only consider unpaired electrons
n_unpaired = 0;
S_elec = R*log(2*n_unpaired+1);
Cp_elec = 0;

%% sum up all contributions
S_total = S_trans+S_rot+S_vib + S_elec;
Cp_total = Cp_trans + Cp_rot + Cp_vib + Cp_elec;

