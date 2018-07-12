function basis =make_global_phi(basis,msh)

% returns the global shape functions and their gradient 
%input ---data, mesh (msh), shape function info (basis)
% Copyright (c) 2018 Ravi Sundaria / Aalto University
Ne= size(msh.t,2);
w=basis.w;
gradPhi_ref=basis.gradPhi_ref;
for k_element = 1:Ne
    [B,~] = get_ElementwiseMapping(msh, k_element);
    basis.B{k_element}=abs(det(B));  
    t=1;
     for i=1:size(w,2)
         gradPhi = (B') \gradPhi_ref(:,t:t+1)';
         basis.gradPhi{i,k_element}=gradPhi;
         t= t+2; 
     end 
end 
basis.gmat=cell2mat(basis.gradPhi);
basis.det=cell2mat(basis.B); 
end 

function [B,b] = get_ElementwiseMapping(msh, el)


B = [msh.p(:,msh.t(2,el))-msh.p(:,msh.t(1,el)) msh.p(:,msh.t(3,el))-msh.p(:,msh.t(1,el))];
b = msh.p(:,msh.t(1,el));


end