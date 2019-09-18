function F = Compute_F( dubar, du, Y, Sigma)
% For linear advection equation:
% Similar to compute_Lu for DO 
% F = -c[ubar_x +  (u_i)_x \Sigma_ij Y_j]

global vmean
global sigma
global xr

N  = size(Y,2);
c = pi*(vmean + sigma*xr);

% terms without Y_i's
F   = -dubar*c';
tmp = zeros(size(du,1), size(c,1));

for i=1:N
    for j=1:N
        tmp = tmp - du(:,i)*Sigma(i,j)*(Y(:,j).*c)';
    end
end
F = F + tmp;
end