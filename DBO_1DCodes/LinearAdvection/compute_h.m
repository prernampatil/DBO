function h = compute_h(Lu, ELu, u, w)
% Compute h(k,j) = h_j(xi_k) = <L[u]-E[L[u]], u_j>(xi_k), k=1,...,Nr, j=1,...,N.

Nr = size(Lu,2);
N  = size(u,2);
h  = zeros(Nr,N);

igd = transpose(Lu - repmat(ELu,1,Nr));		% size(igd) = [Nr Ns]

for j=1:N
    h(:,j) = igd*(u(:,j).*w);
end

end
