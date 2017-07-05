% computation for GeH2
% number of atoms
N = 3;
% pressure
p = 1e5;
% parameters
mu = 74.64;
% parameters
% rotational constants (Hz)
A = 205.0e9;
B = 193.6e9;
C = 99.6e9;
sigma = 2;
% parameters 
% vibrational frequencies
% length of mu = 3*N - 6 fo non-linear molecule, where N is the number of
% atoms
% unit: cm^-1
% after testing, the data difference here does not affect results much
nu = [817 2026 2117]; % data from Wang & Zhang, 2004 
% nu = [935.3 1880.6 1888.0]; % data from Polino et al. 2010
% nu = [932 1910.3 1913.3]; % data again from Polino et al. 2010
% experimental results, consider scale factor = 1
sf_entropy = 1.0;
sf_cp = 1.0;

% number of unpaired electrons
n_unpaired = 0;

% sample T
T = linspace(200,2000,20);
N_sample = length(T);
S = zeros(1,N_sample);
Cp= zeros(1,N_sample);
for i=1:N_sample
    [S(i),Cp(i)] = ...
    thermodynamic_computation(N,T(i),p,mu,A,B,C,sigma,nu,sf_entropy,sf_cp,...
    n_unpaired);
end

%plot(T,S);
%figure
%plot(T,Cp);

% output
fID = fopen('GeH2_tabulated.txt','w');
fprintf(fID,'T (K)\n');
fprintf(fID,'%d ',T);
fprintf(fID,'\n');
fprintf(fID,'Cp (J/mol/K)\n');
fprintf(fID,'%f ',Cp);
fprintf(fID,'\n');
fprintf(fID,'S (J/mol/K)\n');
fprintf(fID,'%f ',S);

