function Sf = sparseFinalize(S, varargin)
% Copyright (c) 2013-2016 Antti Lehikoinen / Aalto University
ri = S.ri-1;
 
Sf = sparse(S.I(1:ri), S.J(1:ri), S.E(1:ri), varargin{:});

end