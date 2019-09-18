function lu  = compute_Lu( dubar, du, Y)
% For the Linear advection equation:
% L[u] = -c[ubar_x + Y_i (u_i)_x]

global vmean
global sigma
global xr
N =size(du,2);

% terms without Y_i's
c = pi*(vmean + sigma*xr);
lu = -dubar*c';
tmp = zeros(size(du,1), size(c,1));
for i=1:N
    tmp = tmp -du(:,i)*(Y(:,i).*c)';
end
lu = lu + tmp;

end

