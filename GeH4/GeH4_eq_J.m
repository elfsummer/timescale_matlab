% load data GeH4 data from CEA runs
% along Jupiter's adiabat
function [X_GeH4_eq] = GeH4_eq_J(T0)
% X_H2 = 0.864, q_H2S = 8.90e-5 X_Ge = 4.4 times solar 
% is hard coded in the data

directory = '../CEA_run/data/Jupiter/new_Ge/4.4_solar/';
fid = fopen(strcat(directory,'Jupiter_1.plt'),'r');
% specify number of t,p+species
n_column = 17;

if fid>=0
   fscanf(fid,'%s',1);
   s = cell(n_column);
   for i=1:n_column
       s{i} = fscanf(fid,'%s',1);
   end
   data_1 = fscanf(fid,'%f');
   data_1 = reshape(data_1,n_column,49);
   data_1 = data_1';
   data_short_1 = data_1(1:11:end,:);
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

T = data_short(:,1);
P = data_short(:,2);
X_GeH4 = data_short(:,10);

X_GeH4_eq = interp1(T,X_GeH4,T0,'cubic');

%legend('Ge', 'GeBr', 'GeBr2', 'GeCl', 'GeCl2'...
%    , 'GeF', 'GeF2', 'GeH4', 'GeO', 'GeS', 'GeS2',...
%    'Ge2', 'Ge(cr)', 'GeO2(II)')