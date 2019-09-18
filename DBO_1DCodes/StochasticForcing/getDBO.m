% MS(:, j) = u(:,t; \xi_j) from PCM solution
% mean = \bar{u}(:,t)
%
% Construct DBO decomposition from the above solution and return
% corresponding u, Sigma and Y for DO decomposition.
%
function [u0, Sigma0, y0] = getDBO(U0, Y0, wr, ~)

[~, N] = size(U0); 

% Find the covariance matrix 
C0 = zeros(N,N);
for i=1:N
    for j=1:N
        C0(i,j) = sum(Y0(:,i) .* Y0(:,j) .* wr);
    end
end

[E,D,~] = svd(C0);

u0=U0;
Sigma0 =E*sqrt(D); 
y0 =inv(sqrt(D))*E'*Y0';

y0= y0';

end
