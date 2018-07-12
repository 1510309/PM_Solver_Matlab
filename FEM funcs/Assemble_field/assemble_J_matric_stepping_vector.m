function [P,S]= assemble_J_matric_stepping_vector(data,msh,Areal,basis,BH_data)

% returns the Jacobian and stiffness matrix 
%input ---data, mesh (msh), shape function info (basis), vector
%potenial(Areal), BH data
% Copyright (c) 2018 Ravi Sundaria / Aalto University
%---------------------------

% define shape functions according to order of elememt
w=basis.w;
[flux_real, ~,flux_abs]=calculate_flux_density2(msh,basis,Areal,[]); 
flux_abs=flux_abs/data.filling;
[reluctivityInElement,reluctivitygradientInElement]= Calculate_reluctivity(data,msh,flux_abs, w,BH_data); 
S = assemble_matrix_ravi2(reluctivityInElement, msh,[],basis);
P = assemble_Jacobian_ravi('grad', 'nodal', 'grad', 'nodal', reluctivityInElement, [], msh,[],basis,reluctivitygradientInElement,flux_real, flux_real);

end
 

function [reluctivityInElement,reluctivitygradientInElement]= Calculate_reluctivity(data,msh,fluxdensity,w,BH)
mu0= (4*pi*1e-7);

% calculate the reluctivity and gradient of reluctivity at each element and

% initialization with vacumme data
reluctivityInElement=(1/mu0)* ones (size(w,2),size(msh.t,2));
reluctivitygradientInElement=zeros (size(w,2),size(msh.t,2));
%for iron materi
ind_iron= find (data.imat==4);
reluctivityInElement(:,ind_iron)=fnval(BH.pp_iron,fluxdensity(:,ind_iron).^2);
reluctivitygradientInElement(:,ind_iron)=fnval(BH.pp_der_iron,fluxdensity(:,ind_iron).^2);
% for shaft material
ind_shaft= find (data.imat==1);
reluctivityInElement(:,ind_shaft)=fnval(BH.pp_shaft,fluxdensity(:,ind_shaft).^2);
reluctivitygradientInElement(:,ind_shaft)=fnval(BH.pp_der_shaft,fluxdensity(:,ind_shaft).^2);
% for PM material
reluctivity_pm=-data.Hc/data.Br;
ind_PM= intersect(find(data.imat>300),find(data.imat<500));
reluctivityInElement(:,ind_PM)= reluctivity_pm;

end