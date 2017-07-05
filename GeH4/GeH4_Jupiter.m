clear all;
%close all;


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

N_K_eddy = 30; 
K_eddy = 10.^linspace(5,12,N_K_eddy);
X_GeH4_quench = zeros(1,N_K_eddy);
for i=1:N_K_eddy
    [quench_T quench_PG] = XGeH4_quench(T_data, PG_data,...
    K_eddy(i), X_H2S, X_H2, g, alpha);
    X_GeH4_quench(i) = GeH4_eq_J(quench_T);
    % X_GeH4_quench(i) = GeH4equ(quench_T,quench_PG,X_H2,X_H2S,X_Ge);
end


% for a variety of K
figure()
set(gca,'FontSize',14);
semilogx(K_eddy,X_GeH4_quench,'b','LineWidth',3);
xlabel('$K_{\rm eddy}$ (cm$^2$ s$^{-1}$)','interpreter','latex')
ylabel('mole fraction','interpreter','latex')
set(gca,'yscale','log');
xlim([1e5,1e12]);

%% add two vertical lines
N_points = 30;
x1 = 1e7*ones(1,N_points);
x2 = 1e9*ones(1,N_points);
y = logspace(-11,-8,N_points);
hold on;
loglog(x1,y,'k--','LineWidth',1);
hold on;
loglog(x2,y,'k--','LineWidth',1);

%% add two horizontal lines

N_points = 30;
x = logspace(4,12,N_points);
y1 = 9e-10*ones(1,N_points);
y2 = 5e-10*ones(1,N_points);
hold on;
loglog(x,y1,'k--','LineWidth',1);
hold on;
loglog(x,y2,'k--','LineWidth',1);

% Create arrow
annotation('arrow',[0.3875 0.45],...
    [0.741857142857143 0.740476190476191]);

% Create arrow
annotation('arrow',[0.591071428571429 0.530357142857143],...
    [0.318047619047619 0.319047619047619]);

% Create arrow
annotation('arrow',[0.307142857142857 0.305357142857142],...
    [0.634714285714286 0.590476190476191]);

% Create arrow
annotation('arrow',[0.703571428571428 0.701785714285714],...
    [0.570428571428572 0.623809523809524]);

% Create textbox
annotation('textbox',...
    [0.214285714285714 0.764285714285715 0.144642857142858 0.0894761904761922],...
    'Interpreter','latex',...
    'String',{'GeH$_4$'},...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation('textbox',...
    [0.714642857142856 0.182380952380954 0.144642857142858 0.0894761904761922],...
    'Interpreter','latex',...
    'String',{'Jupiter'},...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);


%hold on;
%semilogx(K_eddy,X_GeH4_quench,'b--','LineWidth',2);

