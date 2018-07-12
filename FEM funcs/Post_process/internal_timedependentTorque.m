function T = internal_timedependentTorque(A, msh, ri, ro, rotorAngles,basis)
%internal_timedependentTorque torque computation for time-stepping
%problems.
% 
% Copyright (c) 2016 Antti Lehikoinen / Aalto University
% Copyright (c) 2018 Ravi Sundaria / Aalto University
mu0 = pi*4e-7;
msh_comp = struct('p', msh.p, 't', msh.bandData.t_ag);


Nsamples = numel(rotorAngles);

T = zeros(1, Nsamples);

X_quad=basis.IP_cord;% cordinates of integration points
W_quad=basis.w; % weights of integration points
N_quad = numel(W_quad);

for ksample = 1:Nsamples
    if isobject(msh.bandData)
        %[t_local, pnew, t_global] = msh.bandData.t_ag(rotorAngles(ksample));
        t_ag = msh.bandData.t_const; %computing torque from the constant part of ag mesh only
        pnew = msh.bandData.p_virt;
        %msh_comp.t_global = t;
    
    else
        [t_ag, ~, pnew] = updateRotorPosition(rotorAngles(ksample), msh);
    end
    msh_comp.t = t_ag;
    msh_comp.p = pnew;
    
    [F, F0] = mappingTerms(msh_comp);
    detF = mappingDeterminant(F);
    
    for k_quad = 1:N_quad
        B = zeros(2, size(msh_comp.t,2));
        x_global = mappingTimesVector(X_quad(:,k_quad), false, false, F) + F0;
        
        
        %computing flux density
        for kn = 1:size(msh_comp.t,1)
            if isfield(msh, 'symmetrySectors')
                A_true = transpose(A( msh.bandData.el_table(2, msh_comp.t(kn, :)), ksample)) ...
                    .* msh.bandData.el_table(3, msh_comp.t(kn, :));
            else
                A_true = transpose( A(msh_comp.t(kn, :), ksample ) );
            end
            B = B + bsxfun(@times, evaluate_ShapeFunction('curl', 'nodal', kn, F, detF,basis,k_quad), ...
                A_true);
        end
        
        Br_r = dotProduct(B, x_global); %Br * |r|
        Bphi_r = dotProduct(B, [0 -1;1 0]*x_global); %Bphi * |r|

        T(ksample) = T(ksample) + W_quad(k_quad) * sum(( Br_r .* Bphi_r ) ...
            ./ sum(x_global.^2,1).^0.5 .*abs(detF));    
        
    end
end

T = T / (mu0*(ro-ri));

end