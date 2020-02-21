function G = compute_G(p, u, w)
% Compute G(i,j) = <E[L[u]Y_j], u_i> = <p_j, u_i>, i,j=1,...,N
N  = size(u,2);
G  = zeros(N,N);

for i=1:N
    for j=1:N
        G(i,j) = sum(p(:,j).*u(:,i).*w);
    end
end

end
