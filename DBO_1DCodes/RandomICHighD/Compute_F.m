function F = Compute_F(ubar, u, dubar,ddubar, du, ddu, Y, Sigma)
% compute F for the Burgers equation
%	F = -u u_x + nu*u_{xx} + f(x,t;w)
%	  = -ubar ubar_x - Y_j Sigma_ij (u_i ubar)_x - Sigma_ij Sigma_mn Y_j Y_n u_i (u_j)_x +
%		nu*(ubar_{xx} + Y_i (u_i)_{xx}) + f(x,t;w)
%
N  = size(Y,2);
Nr = size(Y,1);
global nu
% terms without Y_i's
F = repmat( -ubar.*dubar + nu*ddubar, 1, Nr );

% for i=1:N
%     for j=1:N
%         % terms with Y_i Y_j's
%         for m=1:N
%             for n=1:N
%                 F = F - (u(:,i).*du(:,m)) * (((Y(:,j).*Y(:,n)))'*Sigma(i,j)*Sigma(m,n));
%             end
%         end
%         % terms with Y_i's
%         tmp = -du(:,i).*ubar - u(:,i).*dubar + nu*ddu(:,i);
%         F = F + tmp*Sigma(i,j)*Y(:,j)';
%     end
% end

A = du*Sigma*Y';
B = u*Sigma*Y';
F = F - A.*B;
temp =  (-du.*ubar - u.*dubar + nu*ddu);
F = F + temp*Sigma*Y';

end