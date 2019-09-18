function [rhs_ubar, rhs_u, rhs_Y, rhs_Sigma ] = compute_rhs_dbo(ubar, u, Y, Sigma, mu, xr, wr, wp, t)
% Compute RHS for mean, Y evolution and Sigma evolution
nu = mu;
[~,N] = size(u);
Ns = size(u,1);

% Calculate the x derivatives
dubar   = fourdifft(ubar,1);
ddubar  = fourdifft(ubar,2);
du      = Diff_rk4(u,1);
ddu     = Diff_rk4(u,2);
C       = Sigma*Sigma';

%  Forcing 
f 	= (0.5*(1+xr)*sin(2*pi*t))';
f	= repmat(f,Ns,1);
% Compute F
F = Compute_F(ubar, u, dubar,ddubar, du, ddu, Y , Sigma, nu, f);

EFu = F*wr;
% MEAN
rhs_ubar = EFu;
F = F - EFu;
% BASIS
rhs_u = zeros(size(u));
for k=1:N
    FYk = F.*Y(:,k)';
    EFYk = FYk*wr;
    rhs_u(:,k) = EFYk;
    for i =1:N
        rhs_u(:,k) = rhs_u(:,k) - sum(EFYk.*u(:,i).*wp)*u(:,i);
    end
end

rhs_u = rhs_u*inv(Sigma); 

% STOCHASTIC BASIS
rhs_Y = zeros(size(Y));
for k =1:N 
    ukF = sum(u(:,k).*F.*wp);
    rhs_Y(:,k) = ukF';
    for i =1:N 
        rhs_Y(:,k) = rhs_Y(:,k) - sum(ukF'.*Y(:,i).*wr)*Y(:,i);
    end
end

rhs_Y = rhs_Y*inv(Sigma');
% SIGMA
rhs_Sigma= zeros(size(Sigma));
for k =1:N
    FYk = F.*Y(:,k)';
    EFYk = FYk*wr;
    for m= 1:N 
        rhs_Sigma(m,k) = sum(EFYk.*u(:,m).*wp);
    end
end

end