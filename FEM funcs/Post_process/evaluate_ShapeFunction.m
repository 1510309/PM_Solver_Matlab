function V = evaluate_ShapeFunction(op, fun, k, varargin)
% 
% Copyright (c) 2016 Antti Lehikoinen / Aalto University
if strcmp(fun, 'nodal')
   phi_ref = {  [-1 -1 1]; [1 0 0]; [0 1 0] };
   if strcmp(op, '')
        V = phi_ref{k}*X;
    elseif strcmp(op, 'grad')
        F = varargin{1}; detF = varargin{2};
        V = mappingTimesVector(phi_ref{k}(1,1:2)', true, true, F, [], detF);
    elseif strcmp(op, 'curl')
        F = varargin{1}; detF = varargin{2};basis=varargin{3};k_quad=varargin{4};
        V = [0 1;-1 0] * mappingTimesVector(basis.gradPhi_ref(k,2*k_quad-1:2*k_quad)', true, true, F, [], detF);
    else
        error('Invalid operator.');
    end
else
    error('Shape function not implemented.');
end


end