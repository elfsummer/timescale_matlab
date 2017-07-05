clear all;
%close all;

%addpath(' /Users/wangdong/Documents/2015_summer/research/PT_profile/Jupiter');

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
%PT_data = csvread('PT_Jupiter.csv');
%T_data = PT_data(:,1);
%PG_data = PT_data(:,2);
X_H2 = 0.864;
enrich_factor = 0.55;
X_As = 4.0e-10*X_H2*enrich_factor;
g = 24.79;
alpha = 1.0;

addpath('./As_equilibrium')
filename = 'As_eq_Jup.dat';
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

% for a variety of K

figure()
set(gca,'FontSize',14);
semilogx(K_eddy,X_AsH3/X_H2,'b','LineWidth',3);
xlabel('$K_{\rm eddy}$ (cm$^2$ s$^{-1}$)','interpreter','latex')
ylabel('mixing ratio','interpreter','latex')
xlim([1e4,1e12])
% add two horizontal lines
hold on;
semilogx(K_eddy,1.1e-10*ones(1,N_K_eddy),'k--','LineWidth',1);
hold on;
semilogx(K_eddy,3.3e-10*ones(1,N_K_eddy),'k--','LineWidth',1);

% add two vertical lines
semilogx(1e7*ones(1,100),linspace(0,1e-9,100),'k--','LineWidth',1);
semilogx(1e9*ones(1,100),linspace(0,1e-9,100),'k--','LineWidth',1);

% Create arrow
annotation('arrow',[0.435714285714286 0.5],...
    [0.663285714285714 0.664285714285714]);

% Create arrow
annotation('arrow',[0.623214285714286 0.558928571428571],...
    [0.439476190476191 0.438095238095238]);

% Create textbox
annotation('textbox',...
    [0.194642857142857 0.739476190476191 0.173214285714286 0.1],...
    'Interpreter','latex',...
    'String',{'AsH$_3$'},...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation('textbox',...
    [0.703571428571428 0.701380952380955 0.132142857142857 0.1],...
    'Interpreter','latex',...
    'String',{'Jupiter'},...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create arrow
annotation('arrow',[0.732142857142857 0.730357142857143],...
    [0.203761904761905 0.252380952380952]);

% Create arrow
annotation('arrow',[0.346428571428571 0.346428571428571],...
    [0.36804761904762 0.326190476190476]);



%hold on;
%semilogx(K_eddy,X_AsH3/X_H2,'b--','LineWidth',2);





