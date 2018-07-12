
% post processing of the simulation 
% Copyright (c) 2018 Ravi Sundaria / Aalto University

% flux density 
figure(4); clf;
drawFluxDensity(msh,A(:,13),wm*13*del_tstp, 'Linestyle', 'none'); colormap('jet'); colorbar;
drawFluxLines(msh,A(:,13), 20,wm*13*del_tstp, 'k');caxis([0 2]); axis equal
title('Flux density');

% torque
% torque computation
figure(5); clf ;
hold on
plot((201:400)*1/50/200,Torque_t(:,201:400),'-*')
xlabel('Time (sec)'); ylabel('Torque (Nm)')
avg_torque= mean(Torque_t(:,201:400));
fprintf('Average Torque %f Nm \n', avg_torque)

% currents phase and line 
figure(6); clf;
hold on
plot((201:400)*1/50/200, data.N_pl*K'*IS(:,201:400),'-*')
xlabel('Time (sec)'); ylabel('Phase currents (A)')
title('Phase Currents');


