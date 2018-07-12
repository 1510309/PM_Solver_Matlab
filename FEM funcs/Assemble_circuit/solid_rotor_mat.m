
function  M= solid_rotor_mat(data, msh, basis)
% function to include the effect of conductivity of  solid part of rotor / 
% magnets 

%input ---data, mesh (msh), shape function info (basis)
% Copyright (c) 2018 Ravi Sundaria / Aalto University


% mass matrix
w=basis.w ;% integration point weight
Phi_ref= basis.Phi_ref;
M= zeros(size(msh.p,2));
Ne= size(msh.t,2);

% initialisation
cond_pm=data.cond_pm;
ind_mat= unique(data.imat);% all mat ind
ind_pm= ind_mat(ind_mat>300);% pm material

for k_element = 1:Ne
    [B,~] = get_ElementwiseMapping(msh, k_element);
    indices = msh.t(:,k_element);
    for j= 1: length(ind_pm)
        if data.imat(k_element)==ind_pm(j)
            for i=1:size(w,2)             
                M_local_vect= Phi_ref(i,:)'*Phi_ref(i,:);
                M(indices, indices) =  M(indices, indices) + ...
                    w(i) * cond_pm* M_local_vect * abs(det(B));
            end   
        end
    end
end

M= sparse(M);

end
function [B,b] = get_ElementwiseMapping(msh, el)
%returns the mapping from the reference element to the global element el,
%such that g_global = B*x_ref + b

B = [msh.p(:,msh.t(2,el))-msh.p(:,msh.t(1,el)) msh.p(:,msh.t(3,el))-msh.p(:,msh.t(1,el))];
b = msh.p(:,msh.t(1,el));

end