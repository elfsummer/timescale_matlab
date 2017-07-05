function [XPH3_quench, quench_T, quench_PG, k_1_f,T] = XPH3_quench(T_data, PG_data, ...
    K_eddy, X_H2O, X_H2, X_P, g, alpha)
% compute quenched abundance of PH3 for a given adiabat

T = T_data(1):0.1:T_data(end);
PG = interp1(T_data,PG_data,T,'cubic');

Rgas = 8.3145e-3;
KBOLTZ=1.38e-16;

% compute t_chem

% rate limiting step assumption:
% PO2 + H2O - H + HOPO2

% nasa polynomial for H
H_nasa_coeff_low = [2.500000000E+00,   7.053328190E-13, ...
               -1.995919640E-15,   2.300816320E-18,  -9.277323320E-22,...
                2.547365990E+04,  -4.466828530E-01];
H_nasa_coeff_high = [2.500000010E+00,  -2.308429730E-11, ...
                1.615619480E-14,  -4.735152350E-18,   4.981973570E-22,...
                2.547365990E+04,  -4.466829140E-01];
G_0_H = evaluate_NASA_poly(H_nasa_coeff_low,H_nasa_coeff_high,T);

% nasa polynomial for PH3
PH3_nasa_coeff_low = [4.17009763E+00, -5.06487157E-03, ...
               2.86027846E-05, -3.13123782E-08, 1.13447768E-11, ...
              2.03144445E+02, 2.02004617E+00  ];
PH3_nasa_coeff_high = [ 3.71229298E+00, 5.85959002E-03, ...
              -2.16607791E-06, 3.56195511E-10, -2.15913467E-14, ...  
              -1.88863997E+02, 1.92781913E+00];
G_0_PH3 = evaluate_NASA_poly(PH3_nasa_coeff_low,PH3_nasa_coeff_high,T);

% nasa polynomial for H2
H2_nasa_coeff_low = [2.344331120E+00,   7.980520750E-03, ... 
               -1.947815100E-05,   2.015720940E-08,  -7.376117610E-12, ...
               -9.179351730E+02,   6.830102380E-01];
H2_nasa_coeff_high = [ 3.337279200E+00,  -4.940247310E-05, ... 
                4.994567780E-07,  -1.795663940E-10,   2.002553760E-14, ...
               -9.501589220E+02,  -3.205023310E+00];
G_0_H2 = evaluate_NASA_poly(H2_nasa_coeff_low,H2_nasa_coeff_high,T);

% nasa polynomial for H2O
H2O_nasa_coeff_low = [4.198640560E+00,  -2.036434100E-03, ...
                6.520402110E-06,  -5.487970620E-09,   1.771978170E-12,...
               -3.029372670E+04,  -8.490322080E-01];
H2O_nasa_coeff_high = [ 3.033992490E+00,   2.176918040E-03, ...
               -1.640725180E-07,  -9.704198700E-11,   1.682009920E-14, ...
               -3.000429710E+04,   4.966770100E+00];
G_0_H2O = evaluate_NASA_poly(H2O_nasa_coeff_low,H2O_nasa_coeff_high,T);

% nasa polynomial for PO2 
PO2_nasa_coeff_low = [3.70338274E+00, 1.99856635E-03, ...
               8.62444075E-06, -1.34479846E-08, 5.58470721E-12, ...
               -3.51050019E+04, 8.53959173E+00];
PO2_nasa_coeff_high = [5.27587261E+00, 1.80638038E-03, ...
              -7.50028350E-07, 1.40062873E-10, -9.17506743E-15, ...
              -3.57348111E+04,-5.74241509E-01 ];
G_0_PO2 = evaluate_NASA_poly(PO2_nasa_coeff_low,PO2_nasa_coeff_high,T);

% nasa polynomial for HOPO2 
HOPO2_nasa_coeff_low = [1.63070941E+00, 2.49990293E-02,  ...
               -2.23292246E-05, 5.83782168E-09, 1.06880806E-12, ...
               -8.72608520E+04, 1.81667382E+01];
HOPO2_nasa_coeff_high = [9.07543004E+00, 3.08252827E-03, ...
                -1.12339846E-06, 1.83697534E-10, -1.11135239E-14, ...
                -8.91864155E+04, -1.97912936E+01 ];
G_0_HOPO2 = evaluate_NASA_poly(HOPO2_nasa_coeff_low,HOPO2_nasa_coeff_high,T);

K_eq_PO2 = exp(-(G_0_PO2 + 3.5*G_0_H2 - G_0_PH3 - 2*G_0_H2O)./(Rgas*T));

n=PG*1e6./(KBOLTZ*T); % unit molecule/cm^3

% need to take a reverse of the reaction rate
A = 6.022e+23;

% for  PO2 + H2O - HOPO2 + H (1)
k_1_r = 3.16e13*exp(-50./(Rgas*T))/A; % cm^3/molec
% need to compute k_1_f

G_0_reaction_1 = G_0_HOPO2 + G_0_H - G_0_H2O - G_0_PO2;

K_eq_reaction_1 = exp(-G_0_reaction_1./(Rgas*T));

k_1_f = K_eq_reaction_1.*k_1_r;

% with k_1_f, we compute the chemistry timescale
t_chem = PG.^(3/2)*X_H2^(7/2)./(k_1_f*X_H2O^3.*K_eq_PO2.*n);



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
XPH3_quench=PH3equ(T(quench_level),PG(quench_level),X_H2,X_H2O,X_P);
quench_T = T(quench_level);
quench_PG = PG(quench_level);

