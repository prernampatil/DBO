%
% MS(:, j) = u(:,t; \xi_j) from PCM solution
% mean = \bar{u}(:,t)
%
% Construct DBO decomposition from the above solution and return
% corresponding u, Sigma and Y for DO decomposition.
%
function [u0, y0] = getBO(U0, Y0, wr, wp)

[Ns, N] = size(U0); 
[Nr, M] = size(Y0);

% Find the covariance matrix 
C0 = zeros(N,N);
for i=1:N
    for j=1:N
        C0(i,j) = sum(Y0(:,i) .* Y0(:,j) .* wr);
    end
end

[E,D,Et] = svd(C0);

u0 = sqrt(D)*E'*U0';
y0 = inv(sqrt(D))*E'*Y0';
u0 = u0';
y0 = y0';

end
