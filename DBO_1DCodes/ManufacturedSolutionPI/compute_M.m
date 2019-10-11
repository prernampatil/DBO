function ret = compute_M(G, ev)
% compute matrix M.
%
%   INPUT:
%       G - matrix G where G_{ij} = < E[L[u]Y_j], u_i >
%       ev - diag(l_1, l_2, ..., l_N) where l_i = <u_i, u_i>
%   OUTPUT:
%       M where M_{ij} = ( G_{ij} + G_{ji} ) / (-ev(i)+ev(j)), i neq j
%               M_{ii} = 0
%

[m n] = size(G);
ret   = zeros(n,n);
ev1   = diag(ev);     % convert into vector

ret = G+G';
for i=1:m
    for j=1:n
        ret(i,j) = ret(i,j)/(-ev1(i)+ev1(j));
    end
end

ret(1:(n+1):end)=0;

end


