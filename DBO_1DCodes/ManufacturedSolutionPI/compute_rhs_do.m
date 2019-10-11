function [rhs_ubar, rhs_u, rhs_Y, u, Y] = compute_rhs_do(ubar, u, Y, N, Ns, Nr, xr, wr, wp, t, mu)
% compute RHS for ubar, u, and Y given all information
global NModes
global nu

C	= ComputeCovBasis(Y,wr);
EYY_inv = inv(C); 

nu = mu;
NModes  = N;
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
rhs_u = (p - u*G)*EYY_inv;

% STOCHASTIC COEFFICIENT
rhs_Y = h;
end
