% equilibrium plots

clear all;
close all;

directory = '/Users/wangdong/Documents/2015_summer/research/CEA_run/data/Jupiter/P/4_times_solar/';
fid = fopen(strcat(directory,'Jupiter_1.plt'),'r');
% specify number of t,p+species
n_column = 20;

if fid>=0
   fscanf(fid,'%s',1);
   s = cell(n_column);
   for i=1:n_column
       s{i} = fscanf(fid,'%s',1);
   end
   data_1 = fscanf(fid,'%f');
   data_1 = reshape(data_1,n_column,49);
   data_1 = data_1';
   data_short_1 = data_1(1:8:end,:);
end
fclose(fid);

fid = fopen(strcat(directory,'Jupiter_2.plt'),'r');
if fid>=0
   fscanf(fid,'%s',1);
   s = cell(n_column);
   for i=1:n_column
       s{i} = fscanf(fid,'%s',1);
   end
   data_2 = fscanf(fid,'%f');
   data_2 = reshape(data_2,n_column,100);
   data_2 = data_2';
   data_short_2 = data_2(1:11:end,:);
end
fclose(fid);

fid = fopen(strcat(directory,'Jupiter_3.plt'),'r');
if fid>=0
   fscanf(fid,'%s',1);
   s = cell(n_column);
   for i=1:n_column
       s{i} = fscanf(fid,'%s',1);
   end
   data_3 = fscanf(fid,'%f');
   data_3 = reshape(data_3,n_column,25);
   data_3 = data_3';
   data_short_3 = data_3(1:6:end,:);
end
fclose(fid);

data_short = [data_short_1;data_short_2;data_short_3];

set(gca,'FontSize',n_column);
T = data_short(:,1);
plot(data_short(:,3),T,'b-','LineWidth',2);
xlim([1e-20,1e-5])
ylim([450,1500])
set(gca,'xscale','log')
%set(gca,'yscale','log')
xlabel('mole fraction','interpreter','latex')
ylabel('T(K)','interpreter','latex')
grid on;

hold on;
for i=4:n_column
    plot(data_short(:,i),T,'b-','LineWidth',2);
end

% Create textbox
annotation('textbox',...
    [0.610156250000001 0.793285714285715 0.0400781250000036 0.048142857142857],...
    'String',{'PH'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation('textbox',...
    [0.488437500000003 0.832857142857143 0.0365624999999969 0.0475714285714287],...
    'String',{'HPO'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation('textbox',...
    [0.422968750000005 0.641428571428571 0.023906249999995 0.0418571428571428],...
    'String',{'PN'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation('textbox',...
    [0.322656250000002 0.509999999999999 0.0281250000000001 0.04],...
    'String',{'PO'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation('textbox',[0.384375 0.75 0.0281250000000001 0.04],...
    'String',{'PO_2'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation('textbox',...
    [0.246875 0.857142857142857 0.0281250000000001 0.04],'String',{'P_2O_3'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation('textbox',...
    [0.278125 0.74 0.0234375000000004 0.0342857142857143],'String',{'PCl'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation('textbox',...
    [0.198125000000001 0.785714285714286 0.0214062499999991 0.04],...
    'String',{'P_3'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation('textbox',...
    [0.19921875 0.578571428571428 0.0281250000000001 0.04],'String',{'P_4'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation('textbox',...
    [0.341875000000004 0.328571428571429 0.0281250000000001 0.04],...
    'String',{'PS'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation('textbox',...
    [0.135625000000001 0.144714285714287 0.058203125 0.0585714285714286],...
    'String',{'POF_3'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation('textbox',...
    [0.134375 0.291428571428572 0.0398437499999997 0.0347142857142861],...
    'String',{'PF_3'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation('textbox',...
    [0.495781250000001 0.441428571428572 0.0281250000000001 0.04],...
    'String',{'P_2'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation('textbox',...
    [0.611406250000006 0.306142857142858 0.049609375 0.0585714285714286],...
    'String',{'PH_2'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation('textbox',...
    [0.792812500000003 0.137571428571428 0.0665624999999967 0.0585714285714286],...
    'String',{'H_3PO_4(L)'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation('textbox',...
    [0.715781250000006 0.506142857142859 0.049609375 0.0585714285714286],...
    'String',{'HOPO'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation('textbox',...
    [0.868906250000003 0.497571428571428 0.0342187499999969 0.0585714285714286],...
    'String',{'PH_3'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation('textbox',...
    [0.7765625 0.773285714285714 0.0656249999999999 0.06],'Interpreter','latex',...
    'String',{'Jupiter'},...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);


% close all;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% for comparison 1 time solar, 10 times solar, and 20 times solar
% 
% %% ten times solar water
% 
% figure()
% 
% set(gca,'FontSize',n_column);
% 
% T = data_short(:,1);
% % HOPO 
% h_HOPO_10 = plot(data_short(:,3),T,'r-','LineWidth',4,'MarkerSize',8);
% hold on;
% % PH3
% h_PH3_10 = plot(data_short(:,7),T,'b-','LineWidth',4,'MarkerSize',8);
% % % H3PO4(L)
% % plot(data_short(:,9),T,'g-','LineWidth',2);
% xlim([1e-12,1e-5])
% ylim([300,1500])
% set(gca,'xscale','log')
% %set(gca,'yscale','log')
% xlabel('mole fraction','interpreter','latex')
% ylabel('T(K)','interpreter','latex')
% grid on;
% 
% %% one time solar water
% %clear all;
% directory = '../../data/P/1_time_water/';
% fid = fopen(strcat(directory,'Saturn_1.plt'),'r');
% n_column=20;
% if fid>=0
%    fscanf(fid,'%s',1);
%    s = cell(n_column);
%    for i=1:n_column
%        s{i} = fscanf(fid,'%s',1);
%    end
%    data_1 = fscanf(fid,'%f');
%    data_1 = reshape(data_1,n_column,100);
%    data_1 = data_1';
%    data_short_1 = data_1(1:11:end,:);
% end
% fclose(fid);
% 
% fid = fopen(strcat(directory,'Saturn_2.plt'),'r');
% if fid>=0
%    fscanf(fid,'%s',1);
%    s = cell(n_column);
%    for i=1:n_column
%        s{i} = fscanf(fid,'%s',1);
%    end
%    data_2 = fscanf(fid,'%f');
%    data_2 = reshape(data_2,n_column,100);
%    data_2 = data_2';
%    data_short_2 = data_2(1:11:end,:);
% end
% fclose(fid);
% 
% fid = fopen(strcat(directory,'Saturn_3.plt'),'r');
% if fid>=0
%    fscanf(fid,'%s',1);
%    s = cell(n_column);
%    for i=1:n_column
%        s{i} = fscanf(fid,'%s',1);
%    end
%    data_3 = fscanf(fid,'%f');
%    data_3 = reshape(data_3,n_column,25);
%    data_3 = data_3';
%    data_short_3 = data_3(1:6:end,:);
% end
% fclose(fid);
% 
% data_short = [data_short_1;data_short_2;data_short_3];
% 
% hold on;
% 
% T = data_short(:,1);
% % HOPO 
% h_HOPO_1 = plot(data_short(:,3),T,'r--','LineWidth',4,'MarkerSize',8);
% % PH3
% h_PH3_1 = plot(data_short(:,7),T,'b--','LineWidth',4,'MarkerSize',8);
% % % H3PO4(L)
% % plot(data_short(:,9),T,'g--','LineWidth',2);
% 
% 
% %% 30 times solar water
% %clear all;
% directory = '../../data/P/30_time_water/';
% fid = fopen(strcat(directory,'Saturn_1.plt'),'r');
% n_column=20;
% if fid>=0
%    fscanf(fid,'%s',1);
%    s = cell(n_column);
%    for i=1:n_column
%        s{i} = fscanf(fid,'%s',1);
%    end
%    data_1 = fscanf(fid,'%f');
%    data_1 = reshape(data_1,n_column,100);
%    data_1 = data_1';
%    data_short_1 = data_1(1:11:end,:);
% end
% fclose(fid);
% 
% fid = fopen(strcat(directory,'Saturn_2.plt'),'r');
% if fid>=0
%    fscanf(fid,'%s',1);
%    s = cell(n_column);
%    for i=1:n_column
%        s{i} = fscanf(fid,'%s',1);
%    end
%    data_2 = fscanf(fid,'%f');
%    data_2 = reshape(data_2,n_column,100);
%    data_2 = data_2';
%    data_short_2 = data_2(1:11:end,:);
% end
% fclose(fid);
% 
% fid = fopen(strcat(directory,'Saturn_3.plt'),'r');
% if fid>=0
%    fscanf(fid,'%s',1);
%    s = cell(n_column);
%    for i=1:n_column
%        s{i} = fscanf(fid,'%s',1);
%    end
%    data_3 = fscanf(fid,'%f');
%    data_3 = reshape(data_3,n_column,25);
%    data_3 = data_3';
%    data_short_3 = data_3(1:6:end,:);
% end
% fclose(fid);
% 
% data_short = [data_short_1;data_short_2;data_short_3];
% 
% % to make plots easier to read, we only plot HOPO PH3 and H3PO4
% hold on;
% T = data_short(:,1);
% % HOPO 
% h_HOPO_30 = plot(data_short(:,3),T,'r:','LineWidth',4,'MarkerSize',8);
% % PH3
% h_PH3_30 = plot(data_short(:,7),T,'b:','LineWidth',4,'MarkerSize',8);
% % % H3PO4(L)
% % plot(data_short(:,9),T,'g.-','LineWidth',2);
% 
% 
% % Create legend
% legend1 = legend([h_PH3_1,h_PH3_10,h_PH3_30],...
%     'H2O: 1 time solar','H2O: 10 times solar','H2O: 30 times solar ');
% set(legend1,...
%     'Position',[0.1298828125 0.657857142857141 0.225390625 0.146071428571429],...
%     'FontWeight','bold',...
%     'FontSize',16);
% 
% % Create textbox
% annotation('textbox',...
%     [0.8171875 0.622857142857143 0.03828125 0.0547142857142856],...
%     'String',{'PH_3'},...
%     'FontWeight','bold',...
%     'FontSize',20,...
%     'FontName','Times New Roman',...
%     'FitBoxToText','off',...
%     'EdgeColor',[1 1 1]);
% 
% % Create textbox
% annotation('textbox',...
%     [0.171250000000001 0.844285714285714 0.03828125 0.0547142857142856],...
%     'String',{'Saturn'},...
%     'FontWeight','bold',...
%     'FontSize',20,...
%     'FontName','Times New Roman',...
%     'FitBoxToText','off',...
%     'EdgeColor',[1 1 1]);
% 
% % Create textbox
% annotation('textbox',...
%     [0.5215625 0.655714285714285 0.03828125 0.0547142857142856],...
%     'String',{'HOPO'},...
%     'FontWeight','bold',...
%     'FontSize',20,...
%     'FontName','Times New Roman',...
%     'FitBoxToText','off',...
%     'EdgeColor',[1 1 1]);
% 
