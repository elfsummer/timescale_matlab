clear all;
%close all;

%addpath('/Users/wangdong/Documents/2015_summer/research/rate_limiting_step/test_rls_PH3');

T_data = [400 450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 ...
    1350 1400 1450 1500];
PG_data = [34.89 52.02 74.41 102.88 138.33 181.70 234.00 296.30 369.74 455.55 555.04 669.65 800.88 ...
    950.36 1119.83 1311.16 1526.32 1767.42 2036.68 2336.48 2669.30...
    3037.79 3444.73];
X_H2 = 0.881;
X_H2S = 3.76e-4*X_H2;
enrich_factor = 10.0;
X_Ge = 8.93e-9*X_H2*enrich_factor;
g = 10.44;
alpha = 0.11;

% equilibrium GeH4
% N = length(PG_data);
% X_GeH4_eq = zeros(1,N);
% for i=1:N
%     X_GeH4_eq(i)=GeH4_eq_S(T_data(i));
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
    [quench_T quench_PG ] = XGeH4_quench(T_data, PG_data,...
    K_eddy(i), X_H2S, X_H2, g, alpha);
    X_GeH4_quench(i) = GeH4_eq_S(quench_T);
end


% for a variety of K
figure()
set(gca,'FontSize',14);
semilogx(K_eddy,X_GeH4_quench,'b','LineWidth',3, 'DisplayName','Saturn, GeH4');
xlabel('$K_{\rm eddy}$ (cm$^2$ s$^{-1}$)','interpreter','latex')
ylabel('mole fraction','interpreter','latex')
set(gca,'yscale','log');
xlim([1e5,1e12]);
ylim([1e-11,1e-8]);

%% add two vertical lines
N_points = 30;
x1 = 1e6*ones(1,N_points);
x2 = 1e8*ones(1,N_points);
y = logspace(-11,-8,N_points);
hold on;
loglog(x1,y,'k--','LineWidth',1);
hold on;
loglog(x2,y,'k--','LineWidth',1);

%% add two horizontal lines

N_points = 30;
x = logspace(4,12,N_points);
y1 = 2e-10*ones(1,N_points);
y2 = 6e-10*ones(1,N_points);
hold on;
loglog(x,y1,'k--','LineWidth',1);
loglog(x,y2,'k--','LineWidth',1);

% Create textbox
annotation('textbox',...
    [0.694642857142857 0.158523809523811 0.196428571428571 0.0904761904761905],...
    'Interpreter','latex',...
    'String',{'Saturn'},...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation('textbox',...
    [0.305357142857142 0.769047619047619 0.121428571428572 0.0751904761904785],...
    'Interpreter','latex',...
    'String',{'GeH$_4$'},...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1],...
    'LineWidth',2);

% Create arrow
annotation('arrow',[0.285714285714285 0.35],...
    [0.371428571428576 0.369047619047624]);

% Create arrow
annotation('arrow',[0.483928571428572 0.410714285714286],...
    [0.200000000000002 0.197619047619048]);

% Create arrow
annotation('arrow',[0.323214285714286 0.325],...
    [0.587095238095239 0.502380952380953]);

% Create arrow
annotation('arrow',[0.662499999999999 0.658928571428571],...
    [0.463285714285715 0.557142857142857]);

%}
%hold on;
%semilogx(K_eddy,X_GeH4_quench,'b--','LineWidth',2);
