function [quench_T, quench_PG] = XAsH3_quench(T_data, PG_data, ...
    K_eddy, X_H2, X_AsH3, g, alpha)
% compute quenched level of As for a given adiabat
% alpha is the coefficient correcting the effective length scale

T = T_data(1):0.1:T_data(end);
PG = interp1(T_data,PG_data,T,'cubic');

Rgas = 8.3145e-3;
KBOLTZ=1.38e-16;

% compute t_chem
n=PG*1e6./(KBOLTZ*T); % unit molecule/cm^3

% rate limiting step: AsH2 + AsH2 - As2H2 + H2
% unit cm^3/s
k_f = 4.63e-11*T.^0.04.*exp(-16.8./T)*10;
% 10 times slower

% AsH3 = AsH2 + H
% H2
a_H2 = [2.344331120E+00,   7.980520750E-03, ...
               -1.947815100E-05,   2.015720940E-08,  -7.376117610E-12,...
               -9.179351730E+02,   6.830102380E-01 ];
b_H2 = [2.932865790E+00,   8.266079670E-04, ...
               -1.464023350E-07,   1.541003590E-11,  -6.888044320E-16, ...
               -8.130655970E+02,  -1.024328870E+00 ];
[Cp_H2,H_H2,S_H2,G_H2] = evaluate_NASA_poly(a_H2,b_H2,T,'NASA7');
a_AsH3 = [9.44635600e-01, 1.50846900e-02, ...
              1.20169600e-06,  -3.39746500e-08, 2.76765600e-11, ...
              7.45916800e+03, 1.83226800e+01  ];
b_AsH3 = [4.17202200e+00, 4.37132300e-03, ...
            2.17757400e-07, -1.18326400e-09, 4.53637400e-13, ...
              6.88291600e+03, 2.80347700e+00];
[Cp_AsH3,H_AsH3,S_Ash3,G_AsH3] = evaluate_NASA_poly(a_AsH3,b_AsH3,T,'NASA7',600);
a_AsH2 = [3.77894500e+00, 1.75923300e-03, ...
             8.07080700e-07, 2.35876300e-09, -3.04352100e-12, ...
             2.00486200e+04, 1.27297400e+00 ];
b_AsH2 = [3.42830700e+00, 3.18114000e-03, ...
              1.46048400e-07, -7.93714500e-10, 1.69441400e-13, ... 
              2.01028200e+04, 2.90470300e+00];
[Cp_AsH2,H_AsH2,S_AsH2,G_AsH2] = evaluate_NASA_poly(a_AsH2,b_AsH2,T,'NASA7',600);
% nasa polynomial for H
H_nasa_coeff_low = [2.500000000E+00,   7.053328190E-13, ...
               -1.995919640E-15,   2.300816320E-18,  -9.277323320E-22,...
                2.547365990E+04,  -4.466828530E-01];
H_nasa_coeff_high = [2.500000010E+00,  -2.308429730E-11, ...
                1.615619480E-14,  -4.735152350E-18,   4.981973570E-22,...
                2.547365990E+04,  -4.466829140E-01];
[Cp_H,S_H,H_H,G_H] = evaluate_NASA_poly(H_nasa_coeff_low,H_nasa_coeff_high,T,'NASA7');

delta_r_G = G_AsH2 + G_H - G_AsH3;
K_eq = exp(-delta_r_G./(Rgas*T));
% equilbrium cosntant for H2 = H + H
K_eq_H = exp(-(2*G_H-G_H2)./(Rgas*T));
% compute chemical timescale
t_chem = X_H2*K_eq_H.*PG./(n.*k_f*X_AsH3.*K_eq.*K_eq);

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

