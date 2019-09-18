function F = Compute_F(ubar, u, dubar,ddubar, du, ddu, Y, Sigma, nu, f)
% compute F for the Burgers equation
%	F = -u u_x + nu*u_{xx} + f(x,t;w)
%	  = -ubar ubar_x - Y_j Sigma_ij (u_i ubar)_x - Sigma_ij Sigma_mn Y_j Y_n u_i (u_j)_x +
%		nu*(ubar_{xx} + Y_i (u_i)_{xx}) + f(x,t;w)
%
N  = size(Y,2);
Nr = size(Y,1);

% terms without Y_i's
F = repmat( -ubar.*dubar + nu*ddubar, 1, Nr );
% terms with Y_i's
for i=1:N
    for j=1:N
        tmp = -du(:,i).*ubar - u(:,i).*dubar + nu*ddu(:,i);
        F = F + tmp*Sigma(i,j)*Y(:,j)';
    end
end

% terms with Y_i Y_j's
for i=1:N
    for m=1:N
        for j=1:N
            for n=1:N
                F = F - (u(:,i).*du(:,m)) * (((Y(:,j).*Y(:,n)))'*Sigma(i,j)*Sigma(m,n));
            end
        end
    end
end

% Forcing terms
F = F + f;
end