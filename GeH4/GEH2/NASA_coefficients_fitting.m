% compute the coeficients of nasa polynomial using known data for
% delta H_f, cp and S as a function of T
% solve the linear regression problem min(Ax-b)'(Ax-b)
% garantee all the quantities are continuous 
% and first derivatives are not ganranteed to be
% continous

% input
% T: an array of temperatures, including the breaking point
% delta_H_f: at 298 K
% S: entropy as a function of T
% Cp: molar heat capacity as a function of T
delta_H_f = 249.5;
T_ref = 298;
T_b = 1500;
% output: 
% coefficients a1,a2,a3,a4,a5,a6,a7 for low temperature
% coefficients
% coefficients b1,b2,b3,b4,b5,b6,b7 for high temperature

% partition data into two parts: >T_b and <T_b
N = length(T);
index_l = find(T<=T_b);
T_l = T(index_l);
Cp_l = Cp(index_l);
S_l = S(index_l);
N_l = length(index_l);

% for T<=T_b
% construct matrix A
% first part
A1 = [ones(N_l,1) T_l' T_l'.^2 T_l'.^3 T_l'.^4 zeros(N_l,1)];
A2 = [log(T_l') T_l' T_l'.^2/2 T_l'.^3/3 T_l'.^4/4 ones(N_l,1)];
% A3 = [1 T_ref/2 T_ref^2/3 T_ref^3/4 T_ref^4/5 1/T_ref 0.];
A = [A1;A2];
% construct the other side
R = 8.3145;
b = [Cp_l'/R; S_l'/R];
x = A\b;
a1 = x(1);a2 = x(2);a3 = x(3);a4 = x(4);a5 = x(5);a7 = x(6);
% compute a6 using the enthalpy relation at 298 K
a6 = (delta_H_f/(R*1e-3*T_ref)-a1-T_ref/2*a2-T_ref^2/3*a3-T_ref^3/4*a4...
    -T_ref^4/5*a5)*T_ref;



% for T >=T_b
index_h = find(T>=T_b);
T_h = T(index_h);
Cp_h = Cp(index_h);
S_h = S(index_h);
N_h = length(index_h);
% construct matrix A
A1 = [ones(N_h,1) T_h' T_h'.^2 T_h'.^3 T_h'.^4 zeros(N_h,1)];
A2 = [log(T_h') T_h' T_h'.^2/2 T_h'.^3/3 T_h'.^4/4 ones(N_h,1)];
% A3 = [1 T_ref/2 T_ref^2/3 T_ref^3/4 T_ref^4/5 1/T_ref 0.];
A = [A1;A2];
% construct the other side
R = 8.3145;
b = [Cp_h'/R; S_h'/R];

x = A\b;
b1 = x(1);b2 = x(2);b3 = x(3);b4 = x(4);b5 = x(5);b7 = x(6);
% compute b6 using the enthalpy relation at 1000K
b6 = (a1-b1)*T_b + (a2-b2)/2*T_b^2 + (a3-b3)/3*T_b^3 + (a4-b4)/4*T_b^4 ...
    +(a5-b5)/5*T_b^5 + a6;

% write into files
fID = fopen('GeH2_poly_data.txt','w');
fprintf(fID,'low temperature part:\n');
fprintf(fID,'%e,%e,%e,%e,%e,%e,%e\n',a1,a2,a3,a4,a5,a6,a7);
fprintf(fID,'high temperature part:\n');
fprintf(fID,'%e,%e,%e,%e,%e,%e,%e',b1,b2,b3,b4,b5,b6,b7);
fclose(fID);
