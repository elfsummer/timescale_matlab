clear all;
%close all;

%addpath(' /Users/wangdong/Documents/2015_summer/research/PT_profile/Saturn');

T_data = [450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 ...
    1350 1400 1450 1500];
PG_data = [26.11888627 37.28801102  51.47093862...
    69.10266133 90.64326221 116.5812512...
     147.4369918  183.7665559  226.1655185...
     275.2726254  331.7730034  396.4009975...
     469.9423185  553.2361601  647.1768949...
     752.7160803  870.8644035  1002.692785...
    1149.33307  1311.980671  1491.89799...
    1690.413976];
%PT_data = csvread('PT_Saturn.csv');
%T_data = PT_data(:,1);
%PG_data = PT_data(:,2);
X_H2 = 0.881;
enrich_factor = 7.5;
X_As = 4.0e-10*X_H2*enrich_factor;
g = 10.44;
alpha = 1.0;

addpath('./As_equilibrium')
filename = 'As_eq_Sat.dat';
delimiterIn = ' ';
headerlinesIn = 2;
As_eq_data = importdata(filename,delimiterIn,headerlinesIn);
data = As_eq_data.data;

N_K_eddy = 50; 
K_eddy = 10.^linspace(4,12,N_K_eddy);
quench_T = zeros(1,N_K_eddy);
quench_PG = zeros(1,N_K_eddy);
X_AsH3 = zeros(1,N_K_eddy);
for i=1:N_K_eddy
    [quench_T(i) quench_PG(i)] = XAsH3_quench(T_data, PG_data,...
    K_eddy(i), X_As, X_H2, g, alpha);
    X_AsH3(i) = exp(interp1(log(data(:,1)/1e5),log(data(:,3)),log(quench_PG(i))));
    
end


figure()
set(gca,'FontSize',14);
semilogx(K_eddy,X_AsH3/X_H2,'b','LineWidth',3);
xlabel('$K_{\rm eddy}$','interpreter','latex')
ylabel('mixing ratio')

xlim([1e4,1e12]);
ylim([0,5e-9]);

% add two horizontal lines
hold on;
semilogx(K_eddy,2e-9*ones(1,N_K_eddy),'k--','LineWidth',1);
hold on;
semilogx(K_eddy,4e-9*ones(1,N_K_eddy),'k--','LineWidth',1);

% add two vertical lines
semilogx(1e6*ones(1,100),linspace(0,5e-9,100),'k--','LineWidth',1);
semilogx(1e8*ones(1,100),linspace(0,5e-9,100),'k--','LineWidth',1);

% Create textbox
annotation('textbox',...
    [0.648214285714286 0.210904761904763 0.196428571428571 0.0904761904761905],...
    'Interpreter','latex',...
    'String',{'Saturn'},...
    'FontWeight','bold',...
    'FontSize',18,...
    'EdgeColor',[1 1 1]);

% Create arrow
annotation('arrow',[0.614642857142857 0.614642857142857],...
    [0.755190476190477 0.703809523809524]);

% Create arrow
annotation('arrow',[0.369642857142857 0.369642857142857],...
    [0.75852380952381 0.707142857142857]);

% Create arrow
annotation('arrow',[0.384285714285714 0.385714285714286],...
    [0.438523809523812 0.495238095238095]);

% Create arrow
annotation('arrow',[0.325714285714285 0.366071428571429],...
    [0.551857142857146 0.552380952380952]);

% Create arrow
annotation('arrow',[0.634642857142856 0.635714285714286],...
    [0.439952380952383 0.5]);

% Create textbox
annotation('textbox',...
    [0.353571428571428 0.556142857142857 0.173214285714286 0.1],...
    'Interpreter','latex',...
    'String',{'AsH$_3$'},...
    'FontWeight','bold',...
    'FontSize',18,...
    'EdgeColor',[1 1 1],...
    'Color',[0 0 1]);

% Create arrow
annotation('arrow',[0.520714285714283 0.464285714285714],...
    [0.550904761904766 0.552380952380952]);


%hold on;
%semilogx(K_eddy,X_AsH3/X_H2,'b--','LineWidth',2);


