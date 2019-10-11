
function [rhs_ubar, rhs_u, rhs_Y, S] = compute_rhs_bo_v1(ubar, u, Y, N, Ns, Nr, xr, wr, wp, t, mu)
% compute RHS for ubar, u, and Y given all information.
% We utilize the fact of the diagonality of ev!
global NModes
%     global Nr
global nu
global x
global RR
global G
global tol_PI

nu = mu;
NModes = N;

dubar   = fourdifft(ubar,1);
ddubar  = fourdifft(ubar,2);
du      = Diff_rk4(u,1);
ddu     = Diff_rk4(u,2);

ev  = ComputeCovBasis(u,wp);

f   = fs(t);
Lu  = compute_Lu(ubar, u, dubar, du, ddubar, ddu, Y, f, mu);
ELu = Lu*wr;				% Ns*1

p	= compute_p(Lu, Y, wr);
h	= compute_h(Lu, ELu, u, wp);
G	= compute_G(p, u, wp);
M 	= compute_M(G, ev);
S 	= compute_S(G, M, ev);

% MEAN
rhs_ubar = ELu;

% BASIS
rhs_u = u*M + p;

% STOCHASTIC COEFFICIENT
rhs_Y = -Y*S'+h;
rhs_Y = rhs_Y*inv(ev);
RR = sqrt(abs(ComputeCovBasis(rhs_Y,wr)));
end

