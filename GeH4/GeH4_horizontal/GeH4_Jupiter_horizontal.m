clear all;
close all;

addpath('/Users/wangdong/Documents/Cornell/2015_summer/research/disequilibrium/GeH4/');

T_data = [303.88 347.92 402.93 ...
    450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 ...
    1350 1400 1450 1500];
PG_data = [7.0 11.0 18.0 ...
    26.11888627 37.28801102  51.47093862...
    69.10266133 90.64326221 116.5812512...
     147.4369918  183.7665559  226.1655185...
     275.2726254  331.7730034  396.4009975...
     469.9423185  553.2361601  647.1768949...
     752.7160803  870.8644035  1002.692785...
    1149.33307  1311.980671  1491.89799...
    1690.413976];
X_H2 = 0.864;
X_H2S = 8.90e-5*X_H2;
enrich_factor = 4.4;
X_Ge = 8.93e-9*X_H2*enrich_factor;
g = 24.79;
alpha = 0.12;
% equilibrium GeH4
% N = length(PG_data);
% X_GeH4_eq = zeros(1,N);
% for i=1:N
%     X_GeH4_eq(i)=GeH4_eq_J(T_data(i));
% end
% figure()
% set(gca,'FontSize',14);
% semilogx(X_GeH4_eq,T_data,'b','LineWidth',4);
% hold on;
% semilogx(X_Ge-X_GeH4_eq,T_data);
% ylabel('$T$','interpreter','latex')
% xlabel('mole fraction')
% grid on;

load K_eddy_750K.mat
N_K_eddy = length(K_eddy_750K);
X_GeH4_quench = zeros(1,N_K_eddy);
for i=1:N_K_eddy
    [quench_T quench_PG] = XGeH4_quench(T_data, PG_data,...
    K_eddy_750K(i), X_H2S, X_H2, g, alpha);
    X_GeH4_quench(i) = GeH4_eq_J(quench_T);
    % X_GeH4_quench(i) = GeH4equ(quench_T,quench_PG,X_H2,X_H2S,X_Ge);
end

% save data
fid = fopen('horizontal_GeH4.txt','w');
fprintf(fid,'%e, %e\n',[phi*180/pi;X_GeH4_quench]);
fclose(fid);

% for a variety of K
figure()
set(gca,'FontSize',14);
plot(phi*180/pi,X_GeH4_quench,'b','LineWidth',3);
xlabel('$\phi$ (degrees)','interpreter','latex')
ylabel('mole fraction','interpreter','latex')
%set(gca,'yscale','log');
xlim([0,90]);
ylim([1e-10,1e-9]);
grid on;

% Create arrow
annotation('arrow',[0.700357142857142 0.700357142857142],...
    [0.829000000000001 0.756190476190478]);

% Create arrow
annotation('arrow',[0.35 0.35],...
    [0.829952380952381 0.757142857142857]);

% Create arrow
annotation('arrow',[0.359285714285714 0.360714285714286],...
    [0.469476190476191 0.547619047619048]);

% Create arrow
annotation('arrow',[0.706071428571426 0.705357142857142],...
    [0.475666666666669 0.552380952380952]);

% Create textbox
annotation('textbox',...
    [0.217857142857143 0.189476190476192 0.316071428571429 0.1],...
    'Interpreter','latex',...
    'String',{'Jupiter GeH$_4$'},...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);




%% add two horizontal lines
N_points = 30;
x = linspace(0,40,N_points);
y1 = 9e-10*ones(1,N_points);
y2 = 5e-10*ones(1,N_points);
hold on;
plot(x,y1,'k--','LineWidth',2);
hold on;
plot(x,y2,'k--','LineWidth',2);

% hold on;
% semilogx(phi*180/pi,X_GeH4_quench,'b--','LineWidth',2);

