
% time stepping calculation of synchronous machines
% Copyright (c) 2018 Ravi Sundaria / Aalto University


% intialization based on voltage/ current source and star/delta winding
% connections
I_abc=Tabc*Iab;% ab to three phase transformation
V_abc=Tabc*vab;% ab to three phase transformation
tstp_per_period= 200; % number of time steps per period
tstp_total = 400; % total number of time steps
del_tstp= 1/fs/tstp_per_period; % step size of time


if conn==1  % conn = 1 (star) if conn = 2 Delta
    K = [1 0 -1; 0 1 -1];
    Q= [0 0 -1;-1 0 -1];
    I_ph = I_abc(1:2); % phase current
else
    K= eye(3); Q= eye(3);
    I_ph= I_abc;
    x_step= [x_step;I_abc(3)];
end
KKT= K*K';
KT=K';

% form the equations for time step
Dsp = [Ds(:,a1) Ds(:,a2)+h*Ds(:,a3)];
Lb= Xb/(2*pi*fs); % end winding inductance
M_sl= solid_rotor_mat(data, msh, basis);
M_slp = make_periodic(M_sl,a1,a2,a3,h);

% initialize the variables
A= zeros(size(msh.p,2),tstp_total);
IS= zeros(length(I_ph),tstp_total);
maxiter=20;
plp=data.poles/2; % number pole of pole pairs
w=2*pi*fs; % elec speed
wm=w/plp;% mech speed
for tstp = 2: tstp_total
    t= (tstp-1)*del_tstp;
    t1= (tstp-2)*del_tstp; % at previous time step
    v= sqrt(2)*Vs*[cos(w*t) cos(w*t-2*pi/3) cos(w*t-4*pi/3)];
    vt1=  sqrt(2)*Vs*[cos(w*t1) cos(w*t1-2*pi/3) cos(w*t1-4*pi/3)];
    % t1 means t-1 time step
    sys_t1= [-(M_slp)/del_tstp zeros(size(Dsp'*KT));
        -K*Dsp  Lb*KKT/Ns/Lc];
    
    fv= [-J_pm_p; Q*(v')*del_tstp/Ns/Lc];
    
    f_t1= sys_t1*x_step+fv;
    S_ag =get_MovingBandMatrix_ravi_ordr((rot_anglm+wm*t),IP1,basis,msh);% airgap matrix
    
    % calculate the source vector for the time stepping
    
    % the matrix equations
    An= zeros(size(msh.p,2),1);
    % nonlinear iterations
    for iter = 1: maxiter
        An(a1)= x_step(1:length(a1)); % vector potential at this nonlinear step
        An(a2)= x_step(length(a1)+[1:length(a2)]);
        An(a3)=h*An(a2);
        
        [P,St2]= assemble_J_matric_stepping_vector(data,msh,An,basis,BH);
        
        P= P+S_ag;
        St = St2+S_ag;
        Stp = make_periodic(St,a1,a2,a3,h);
        Pp= make_periodic(P,a1,a2,a3,h);
        J_step= [(Pp+M_slp/del_tstp) Dsp'*KT;
            K*Dsp  -(Rs*del_tstp+Lb)*KKT/Ns/Lc];
        
        % step and nth nonlinear itration
        sys_n= [(Stp+M_slp/del_tstp) Dsp'*KT;
            K*Dsp  -(Rs*del_tstp+Lb)*KKT/Ns/Lc];% system matrix of nonlinear system see St
        
        f_t = sys_n*x_step;
        f_step = f_t+f_t1; % f_t depends on nonlin iteration; f_t1 on previous time step value
        del_x= -J_step\f_step;
        x_step= x_step+del_x;
        
        resid=norm(del_x)/norm(x_step);
        %   check convergence with respect to a, ur, Is
        flag =0;
        if resid< convergence_limit %&& norm(del_x(2*length(ar)+(1:2*length(ur_R))))/norm(x(2*length(ar)+(1:2*length(ur_R))))< convergence_limit && norm(del_x(end-7:end))/norm(x(end-7:end))< convergence_limit %,iter>max_iter)
            flag=1;
            break
        end
        if(flag==1)
            break
        end
        
        X= ['timestep=',num2str(tstp),'   ','iteration=',num2str(iter),'   ',  'residual=',num2str(resid)];
        disp(X)
    end
    A(a1,tstp)= x_step(1:length(a1));
    A(a2,tstp)= x_step(length(a1)+[1:length(a2)]);
    A(a3,tstp)= - A(a2,tstp);% magnetic vector potential
    IS(:,tstp)= x_step(end-length(I_ph)+1:end); % stator current
    
    Torque_t(tstp) = internal_timedependentTorque(A(:,tstp), msh,data.ro_out/2,data.st_in/2, (rot_anglm+wm*t),basis);
    
end
Torque_t=Torque_t*msh.symmetrySectors*data.Lc;

