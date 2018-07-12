function J= assemble_Jacobian_ravi(op1, f1, op2, f2, v, elements, msh,S,basis,dnu,flux_real, flux_imag)
% J = Jacobian matrix--> derivative of stiffness with magnetic vector
% potential

% So far only op 'grad' and 'nodal' works for higher order elements
%   S = assemble_matrix(op1, f1, op2, f2, v, elements, msh, S) assembles
%   the sparse matrix struct S with the entries
%   S_ij = Int( v (op1 f1).(op2 f2) ),
%   integrated over the "elements".
%
%   Possible input:
%   op1,op2 = 'grad'
%   f1,f2 = 'nodal'
%   v = either a separate value for each element, or a constant to be used
%   with all elements
%   elements = either a list, or [] to evalute over the entire mesh
%   msh = the mesh struct
%   S = existing matrix struct, or [] to create new ones

%input ---data, mesh (msh), shape function info (basis)
% Copyright (c) 2018 Ravi Sundaria / Aalto University

%no elements listed --> going over all
if ~any(elements)
    elements = 1:size(msh.t,2);
end

%size of v(x) too small --> assuming constant value
if numel(v) == 1
    v = v(1) * ones(1, size(msh.t,2));
end

F = mappingTerms(msh, elements);
detF = mappingDeterminant(F);
DETF = abs(detF);

%numbers of something
Np = size(msh.p, 2);
Ne = numel(elements);
NO_f1 = size(msh.t,1);
NO_f2 = size(msh.t,1);


%defining function handle for (op1 f1)

if strcmp(f1, 'nodal')
    
    INDM_1 = msh.t(:, elements);
    
else
    error(['Undefined function type ' f1]);
end



if strcmp(f2, 'nodal')
    
    INDM_2 = msh.t(:, elements);
else
    error(['Undefined function type ' f1]);
end

%initializing quadrature calculation
E = cell(NO_f1, NO_f2); [E{:}] = deal( zeros(1, numel(elements)) );
N_quad = numel(basis.w);
W_quad=basis.w;
%calculating entries


Br_x=flux_real(1:2:end,:); % the x component
Br_y=flux_real(2:2:end,:); % the y component

Bi_x=flux_imag(1:2:end,:); % the x component
Bi_y=flux_imag(2:2:end,:); % the y component
Br_xy=Br_x.*Bi_y;
Br_xx=Br_x.*Bi_x;
Br_yy=Br_y.*Bi_y;
Br_yx=Br_y.*Bi_x;



if numel(v)==0
    
    Br_xy=2*dnu.*Br_xy;
    Br_xx=2*dnu.*Br_xx;
    Br_yy=2*dnu.*Br_yy;
    Br_yx=2*dnu.*Br_yx;
else 
    
    Br_xy=2*dnu.*Br_xy;
    Br_xx=v+2*dnu.*Br_xx;
    Br_yy=v+2*dnu.*Br_yy;
    Br_yx=2*dnu.*Br_yx;
    
end
%%%%%

for k_quad = 1:N_quad
    for ind_1 = 1:NO_f1
        for ind_2 = 1:NO_f2
        
            
            
            
            V= basis.gmat(2*k_quad-1:2*k_quad,ind_2:NO_f2:end);
            v_curl=[V(2,:);-V(1,:)];
            abcd= [Br_xx(k_quad,:).*v_curl(1,:)+Br_xy(k_quad,:).*v_curl(2,:);
                Br_yx(k_quad,:).*v_curl(1,:)++Br_yy(k_quad,:).*v_curl(2,:); ];
            V2= basis.gmat(2*k_quad-1:2*k_quad,ind_1:NO_f1:end);
            v2_curl=[V2(2,:);-V2(1,:)];
            
            E{ind_1, ind_2} = E{ind_1, ind_2} + W_quad(k_quad) * ...
                dotProduct(  v2_curl, abcd);
        end
    end
end

%adding entries to matrix struct
for ind_1 = 1:NO_f1
    for ind_2 = 1:NO_f2
        S = sparseAdd( INDM_1(ind_1,:), INDM_2(ind_2,:), ...
            E{ind_1,ind_2} .* DETF, S);
    end
end
J = sparseFinalize(S, Np, Np);
end