
% Static simulation of synchronous machine 
% Copyright (c) 2018 Ravi Sundaria / Aalto University

rot_anglm= rot_angl*pi/180/2;% mech rot angle
N_pl=data.N_pl; % number of parallel paths
Xb= data.Xb*N_pl*(fs/50); % satator end winding reactance of a phase at 50 Hz
Rs=  data.Rs*N_pl; % Resistance of stator phase winding
Ns= msh.symmetrySectors/N_pl; % number of symmetry sectors = Number of stator slots/ number of stator slots in solution sector
Lc=data.Lc; % effective length of stator core
Rs=Rs*1.2;% rougly effect of temprature
h= msh.periodicityCoeff; 
S= assemble_S_arrays(data,msh,basis); % Linear material stifness matrix
Ds= assemble_stator_flux(data,msh,basis);
J_pm= assemble_S_pm(data,msh,basis);
%airgap matrix ;% effect of rotation angle added
S_ag =get_MovingBandMatrix_ravi_ordr(rot_anglm,IP1,basis,msh);


vs= sqrt(2)*[Vs Vs*cosd(-120) Vs*cosd(-240)]';% instentaneaous voltage
% park transformation
Tab= (2/3)*[1 exp(1j*2*pi/3) exp(1j*4*pi/3)];
Tab=[real(Tab);imag(Tab)];
Tabc= (3/2)*[2/3 0; -1/3 1/sqrt(3); -1/3 -1/sqrt(3)];

if conn == 1
    vs= vs/sqrt(3);% star connection per phase voltage
    vab= Tab*vs; % conversion of voltage
    vab= [cos(pi/6) sin(pi/6); -sin(pi/6) cos(pi/6)]*vab;% due to line and phase
else
    vab= Tab*vs; % conversion of voltage
end

Ds_ab= Tab*Ds; 
S=S+S_ag; % addting the airgap contribution


% periodic boundaries
Sp= make_periodic(S,a1,a2,a3,h);
Ds_abp= [Ds_ab(:,a1) Ds_ab(:,a2)+h*Ds_ab(:,a3)];
J_pm_p= [J_pm(a1); J_pm(a2)+h*J_pm(a3)];
% circuit equation
iab= [1 0 ; 0 1]; % identity matrix
Pab= [0 -1; 1 0]; % effect of imag part of alpha beta


%----------------------------------------------------------------------
% Linear simulation
sys_mat= [Sp (3/2)*Ds_abp'
    Pab*Ds_abp -(Rs*iab+Xb*Pab)./(  Ns*Lc*(2*pi*fs))  ];
f_mat= [J_pm_p; -vab./(  Ns*Lc*(2*pi*fs))];% source mat
sol= sys_mat\f_mat;

A_pot= zeros(size(S,1),1); % initialize the magnetic vector potential
A_pot(a1)= sol(1:length(a1)); A_pot(a2)= sol(length(a1)+1:end-2); A_pot(a3)= -A_pot(a2);
x_step= sol; % initialization of nonlinear simulation with linear sol
An= zeros(size(A_pot));
%------------------------------------------------------------------------
% Nonlinear simulation
max_iter= 30;
convergence_limit=1e-6;
for iter= 1:20
    
    An(a1)= x_step(1:length(a1)); % vector potential at this nonlinear step
    An(a2)= x_step(length(a1)+[1:length(a2)]);
    An(a3)=h*An(a2);
    % Jacobian matrix
    [P,St2]= assemble_J_matric_stepping_vector(data,msh,An,basis,BH);
    P=P+S_ag;%  adding the contribution of airgap and mass matrix;
    St=St2+S_ag;
    Stp= make_periodic(St,a1,a2,a3,h);
    Pp= make_periodic(P,a1,a2,a3,h);
    
    
    J_step= [Pp (3/2)*Ds_abp'
        Pab*Ds_abp -(Rs*iab+Xb*Pab)./(  Ns*Lc*(2*pi*fs))  ];
    sys_n= [Stp (3/2)*Ds_abp'
        Pab*Ds_abp -(Rs*iab+Xb*Pab)./(  Ns*Lc*(2*pi*fs))  ];
    f_t = sys_n*x_step;% nonlinear iteration f(iter-1)
    f_step = f_t-f_mat;
    del_x= -J_step\f_step;
    x_step= x_step+del_x;
    
    resid=norm(del_x)/norm(x_step); % residual
  
    flag =0;
    if resid< convergence_limit 
        flag=1;
        break
    end
    if(flag==1)
        break
    end
    
    X= ['iteration=',num2str(iter),'   ',  'residual=',num2str(resid)];
    disp(X)
end
An(a1)= x_step(1:length(a1)); % vector potential at this nonlinear step
An(a2)= x_step(length(a1)+[1:length(a2)]);
An(a3)=h*An(a2);

Iab=x_step(end-1:end);
% current

% Terminal quantities
U = norm(vab)/sqrt(2);% rms from peak
I = norm(Iab)/sqrt(2);
theta = atan2(vab(2),vab(1))*180/pi;
if conn == 1
    U = U*sqrt(3); % Star
    theta = theta+30;
else
    I = I*sqrt(3); % Delta
end
Line_I= I*N_pl; %Line current and  effect of parallel path






