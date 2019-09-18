function ret = ComputeCovBasis(u, w)
% compute <u_i, u_j> for all i and j.
%
% INPUT:
% 	u - basis of size [Ns,N]
% 	w - weights of size Ns
% global dY
	[Ns N] = size(u);
	ret 	= zeros(N,N);

	for i=1:N
		for j=1:N
			ret(i,j) = sum( u(:,i) .* u(:,j) .* w);
		end
	end

% 	ret = ret+transpose(tril(ret,-1));	% symmetrize
end
