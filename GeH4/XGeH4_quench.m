function [quench_T, quench_PG] = XGeH4_quench(T_data, PG_data, ...
    K_eddy, X_H2S, X_H2, g, alpha)
% compute quenched level of GeH4 for a given adiabat
% alpha is the coefficient correcting the effective length scale

T = T_data(1):0.1:T_data(end);
PG = interp1(T_data,PG_data,T,'cubic');

Rgas = 8.3145e-3;
KBOLTZ=1.38e-16;

% compute t_chem
n=PG*1e6./(KBOLTZ*T); % unit molecule/cm^3

% rate limiting step: GeH2 + H2S - H2Ge=S + H2
% Fegley and Lodders (1994)
% unit cm^3/s
% uncetainty +-1000 K
% k_f = 1e-11*exp(-6000./T);

% rate limiting step: GeH2 + H2S - HGe=SH + H2
% k_f = 3.57e-14*T.^0.7.*exp(-4956./T);
% faster by factor of 5
k_f = 3.57e-14*T.^0.7.*exp(-4956./T);

% GeH4 = GeH2 + H2
a_GeH2 = [3.776667e+00,9.625977e-04,3.939253e-06,-3.805175e-09,...
    1.022434e-12,2.881190e+04,4.605981e+00 ];
b_GeH2 = [3.776667e+00,9.625977e-04,3.939253e-06,-3.805175e-09,...
    1.022434e-12,2.881190e+04,4.605981e+00];
[Cp_GeH2,H_GeH2,S_GeH2,G_GeH2] = evaluate_NASA_poly(a_GeH2,b_GeH2,T,'NASA7');

% H2
a_H2 = [2.344331120E+00,   7.980520750E-03, ...
               -1.947815100E-05,   2.015720940E-08,  -7.376117610E-12,...
               -9.179351730E+02,   6.830102380E-01 ];
b_H2 = [2.932865790E+00,   8.266079670E-04, ...
               -1.464023350E-07,   1.541003590E-11,  -6.888044320E-16, ...
               -8.130655970E+02,  -1.024328870E+00 ];
[Cp_H2,H_H2,S_H2,G_H2] = evaluate_NASA_poly(a_H2,b_H2,T,'NASA7');

% GeH4
a_GeH4 = [2.54992789E+00 7.13885765E-03 1.43758337E-05   ...
-2.33592977E-08 9.65676013E-12 9.69756465E+03 9.02678812E+00 ];
b_GeH4 = [ 5.41474159E+00 7.24155154E-03 -2.71818301E-06 4.51535021E-10...
    -2.75635275E-14 8.46356611E+03 -7.83419271E+00 ];
[Cp_GeH4,H_GeH4,S_GeH4,G_GeH4] = evaluate_NASA_poly(a_GeH4,b_GeH4,T,'NASA7');

delta_r_G = G_GeH2 + G_H2 - G_GeH4;
K_eq = exp(-delta_r_G./(Rgas*T));

t_chem = (PG*X_H2)./(k_f.*K_eq.*n*X_H2S);

% compute mixing timescale
% calculate pressure scale height
mu = X_H2*2 + (1-X_H2)*4;
% scale height (~km)
H = Rgas*T/(mu*1e-3*g);

L_eff = alpha*H;

t_mix = L_eff.^2*1e10/K_eddy;

%semilogx(t_chem,T,t_mix,T);

%find the quench level
[min_value,quench_level]=min(abs(t_chem-t_mix));
%XGeH4_quench=GeH4equ(T(quench_level),PG(quench_level),X_H2,X_H2S,X_Ge);
quench_T = T(quench_level);
quench_PG = PG(quench_level);

