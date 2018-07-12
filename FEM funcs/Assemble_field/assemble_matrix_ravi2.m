function S_ag = assemble_matrix_ravi2(v, msh,S,basis)
% S = assemble_matrix returns the matrix struct S

%input ---reluctivity (v), mesh (msh), shape function info (basis)
% Copyright (c) 2018 Ravi Sundaria / Aalto University

    elements = 1:size(msh.t,2);


%size of v(x) too small --> assuming constant value
if numel(v) == 1
    v = v(1) * ones(1, numel(elements));
end

% F = mappingTerms(msh, elements);
% detF = mappingDeterminant(F);
% DETF = abs(detF);

%numbers of something
Np = size(msh.p, 2);
NO_f1 = size(msh.t,1);
NO_f2 = size(msh.t,1);






%initializing quadrature calculation
E = cell(NO_f1, NO_f2); [E{:}] = deal( zeros(1, numel(elements)) );
N_quad = numel(basis.w); 
 W_quad=basis.w; 
%calculating entries
for k_quad = 1:N_quad
    for ind_1 = 1:NO_f1
        for ind_2 = 1:NO_f2
            E{ind_1, ind_2} = E{ind_1, ind_2} + W_quad(k_quad) * ...
                v(k_quad,:).*dotProduct( basis.gmat(2*k_quad-1:2*k_quad,ind_1:NO_f2:end), basis.gmat(2*k_quad-1:2*k_quad,ind_2:NO_f2:end) );
        end
    end
end

%adding entries to matrix struct
for ind_1 = 1:NO_f1
    for ind_2 = 1:NO_f2
        S = sparseAdd( msh.t(ind_1,:), msh.t(ind_2,:), ...
             E{ind_1,ind_2} .* basis.det, S);
    end
end
        S_ag = sparseFinalize(S, Np, Np);    
end