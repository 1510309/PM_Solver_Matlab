function S= assemble_S_arrays(data,msh,basis)

% returns the stiffness matrix for linear materials
%input ---data, mesh (msh), shape function info (basis)
% Copyright (c) 2018 Ravi Sundaria / Aalto University

% take the material data of stator and rotor cores from data
% Linear element assumption
mu0= 4*pi*1e-7;
mur= 1000; % relative permeability

%shape functions 
w=basis.w ;% integration point weight
gradPhi_ref=basis.gradPhi_ref;
Ne= size(msh.t,2);

% initialisation
S= zeros(size(msh.p,2));
for k_element = 1:Ne
    [B,~] = get_ElementwiseMapping(msh, k_element);
    indices = msh.t(:,k_element);
    if data.imat(k_element)==4% iron material
        reluctivityInElement=1/(mu0*mur);
    elseif floor(data.imat(k_element)/100)==3% PM material
 
       reluctivity_pm=-data.Hc/data.Br;
       reluctivityInElement=reluctivity_pm;
    else
        reluctivityInElement=1/(mu0);
    end
    t=1;
    for i=1:size(w,2)
        gradPhi = (B') \gradPhi_ref(:,t:t+1)'; 
        S_local_vect =  gradPhi' * gradPhi;
        
        S(indices, indices) =  S(indices, indices) + ...
            w(i) * reluctivityInElement * S_local_vect * abs(det(B));
          t= t+2;
    end
end
end

function [B,b] = get_ElementwiseMapping(msh, el)

B = [msh.p(:,msh.t(2,el))-msh.p(:,msh.t(1,el)) msh.p(:,msh.t(3,el))-msh.p(:,msh.t(1,el))];
b = msh.p(:,msh.t(1,el));

end