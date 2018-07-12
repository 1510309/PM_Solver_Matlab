
% Starter file of simulating permanent magnet machine 

% Copyright (c) 2018 Ravi Sundaria / Aalto University

clc; clear; close all; 
addpath(genpath('data_files'))
addpath(genpath('FEM funcs'))
%input 
% connection ; supply voltage (rms); rotation angle 
rot_angl=55;%input('rotation angle in elec degrees \n'); % rotational angle 
Vs= 400;%input('supply rms voltage \n');
conn=1;%input('connection for stator winding Star(1) Delta(2) \n'); 
fs= 50; % supply frequency 
%---------------------------------------------------------------------------
% Machine data 
% Mesh 
load('Mesh_PM')
% msh field {p= mesh points; t= element; rotel= rotor elements; t_ag= airgap elements}
load('data_machine.mat'); 
load('BH_data'); % magnetization data of iron and shaft 
Machine_data; 
% nodes a2 and a3 are periodic; ab and a1 are dirichlet and nondirichlet nodes
a1=data.a1;a2=data.a2; a3=data.a3; ab=data.ab; 
% Shape functions/ basis functions
node_num= size(msh.t,1);
IP1= 3;%polynomial order for selecting the number of integraion points per element
[Phi_ref, gradPhi_ref, Wt, IP_cord]= calculate_Phiref_gradPhi1 (node_num,IP1);
basis= struct('Phi_ref', Phi_ref, 'gradPhi_ref',gradPhi_ref,'w', Wt,'IP_cord',IP_cord);
basis=make_global_phi(basis,msh);

% For rotation 
msh.bandData = initializeBandData_ordr(msh, [], msh.rotel, msh.t_ag);

% Dc simulation 
PM_DC
% Time-Stepping simulation 
PM_Time_stepping
% post_Processing 
post_process