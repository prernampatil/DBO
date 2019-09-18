function ret = compute_S(G, M, ev)
% compute matrix  = G + M \Lambda
%
% INPUT:
%       G - matrix G where G_{ij} = < E[L[u]Y_j], u_i >
%       M - matrx M that is computed from compute_M function
%       ev - diag(l_1, l_2, ..., l_N) where l_i = <u_i, u_i>
%   OUTPUT:
%       S where S_{ij} = G_{ij} + ev(i)*M_{ij}, i neq j
%               S_{ii} = G_{ii}
%   
    
    n = size(G,1);
    ret = G + ev*M;
    ret(1:(n+1):end)=G(1:(n+1):end);
end
