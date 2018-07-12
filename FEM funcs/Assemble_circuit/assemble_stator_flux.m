function Ds= assemble_stator_flux(data,msh,basis )
%returns stator flux linkage matrix 
%input ---data, mesh (msh), shape function info (basis)
% Copyright (c) 2018 Ravi Sundaria / Aalto University

m=3; % number of phases
Nn= size(msh.p,2);% number of nodes

% shape function info
w=basis.w;
Phi_ref= basis.Phi_ref;

Sn = data.Sn;% ALA1 the area of conductor reigon in a stator slot
Ncn= data.Ncn; % number of conductors in a stator slot

coil_mat_index=[41 42 43;44 45 46];
coil_sign= [1 1 1;-1 -1 -1]; % hard coded coil_mat_index are 41 42 46(46 is for - side of third phase)

% Ds dimension m (phases)* Nn (number of nodes )
Ds = zeros(m,Nn); % m number of phases and Nn are number of nodes in the studied mesh
for ph =1:m
    for sgn = 1:2
        stator_coil_ele= find(data.imat==coil_mat_index(sgn,ph)); % elements belong to phase ph
        for num= 1:length(stator_coil_ele)
            ele= stator_coil_ele(num);
            Ar = get_ElementwiseMapping(msh, ele); % calculate Ar
            indices= msh.t(:,ele);
            Ds(ph,indices)= Ds(ph,indices)-((Ncn/Sn)*coil_sign(sgn,ph))* sum((repmat(w',1,size(Phi_ref,2)).*Phi_ref),1) * abs(det(Ar));
        end
    end
end
end

function [Ar] = get_ElementwiseMapping(msh, el)

Ar = [msh.p(:,msh.t(2,el))-msh.p(:,msh.t(1,el)) msh.p(:,msh.t(3,el))-msh.p(:,msh.t(1,el))];
end