function [rhs_ubar rhs_u rhs_Y u Y] = compute_rhs_do(ubar, u, Y, N, Ns, Nr, xr, wr, wp, t, mu)
% compute RHS for ubar, u, and Y given all information
global NModes
%     global Nr
global nu
global x
global Stoch_forcing
global tol_PI
C	= ComputeCovBasis(Y,wr);
EYY_inv = inv(C); 

nu = mu;
NModes = N;
dubar   = fourdifft(ubar,1);
ddubar  = fourdifft(ubar,2);
du      = Diff_rk4(u,1);
ddu     = Diff_rk4(u,2);

f   = fs(t);
Ef 	= f*wr;					% E[forcing] of size (Ns*1) vector
Lu  = compute_Lu(ubar, u, dubar, du, ddubar, ddu, Y, f, mu);
ELu = Lu*wr;				% Ns*1

p	= compute_p(Lu, Y, wr);
h	= compute_h(Lu, ELu, u, wp);
G	= compute_G(p, u, wp);

% MEAN
rhs_ubar = ELu;

% BASIS
LU = nu*ddu  - repmat(dubar,1,NModes).*u - du.*repmat(ubar,1,NModes);
BU = [];
YY = [];
EfY = [];
for j=1:NModes
    for k=1:NModes
        BU = [BU u(:,j).*du(:,k)];
        YY = [YY Y(:,j).*Y(:,k)];
    end
    fy = f.*repmat(Y(:,j)',Ns,1);
    EfY = [EfY fy*wr];
end

EYYY = Compute_EYYY(Y);
f_basis = [];
for i=1:NModes
    eyyy    =  EYYY(:,:,i);
    f_basis = [f_basis BU*eyyy(:)];
end
quadratic = diag(ComputeCovBasis(f_basis,wp));
Stoch_forcing = diag(ComputeCovBasis(EfY,wp));
f_basis = (EYY_inv*(-f_basis'+EfY'))';
f_basis = LU + f_basis ;
rhs_u = f_basis - u*(u'*(f_basis.*repmat(wp,1,NModes)));

rhs_u = (p - u*G)*EYY_inv;

% STOCHASTIC COEFFICIENT
rhs_Y = h;
end
