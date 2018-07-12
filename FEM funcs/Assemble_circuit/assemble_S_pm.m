function [J_pm]= assemble_S_pm(data,msh,basis)

% returns--- PM source matrix (J_pm)
%input ---data, mesh (msh), shape function info (basis)
% Copyright (c) 2018 Ravi Sundaria / Aalto University

%shape functions and numerical integration for  first order only
w=basis.w ;% integration point weight
gradPhi_ref=basis.gradPhi_ref;

Ne= size(msh.t,2); % number of elements on mesh
Np = size(msh.p,2); % number of nodes on mesh
J_pm = zeros (Np,1);
% initialisation

for k_element = 1:Ne
    [B,~] = get_ElementwiseMapping(msh, k_element);
    indices = msh.t(:,k_element);
    if floor(data.imat(k_element)/100)==3% PM material
        
        midpnt= mean(msh.p(:,msh.t(:,k_element)),2); % centroid of element
        theta= atan(midpnt(2)/midpnt(1));
        
        
        t=1;
        for i=1:size(w,2)
            gradPhi = (B') \gradPhi_ref(:,t:t+1)'; %gradients of shape functions of the GLOBAL element
            
            
            %calculating the contribution with gaussian quadrature
            J_local_vect=gradPhi(1,:)*data.Hc*sin(theta)-gradPhi(2,:)*data.Hc*cos(theta);
            
            J_pm(indices,1)= J_pm(indices,1)+...
                w(i) * J_local_vect' * abs(det(B));
            t= t+2;
        end
    end
end
end

function [B,b] = get_ElementwiseMapping(msh, el)
%returns the mapping from the reference element to the global element el,
%such that g_global = B*x_ref + b
B = [msh.p(:,msh.t(2,el))-msh.p(:,msh.t(1,el)) msh.p(:,msh.t(3,el))-msh.p(:,msh.t(1,el))];
b = msh.p(:,msh.t(1,el));

end