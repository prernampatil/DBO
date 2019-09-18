function [rhs_ubar, rhs_u, rhs_Y] = compute_rhs_do(ubar, u, Y, N, Ns, Nr, xr, wr, wp, t)
% compute RHS for ubar, u, and Y given all information
global x
global vmean
global sigma
C	= ComputeCovBasis(Y,wr);

dubar   = fourdifft(ubar,1);
du = Diff_rk4(u,1);
Lu  = compute_Lu( dubar,du, Y);

ELu = Lu*wr;    % Ns*1

h	= compute_h(Lu, ELu, u, wp); 
p	= compute_p(Lu, Y, wr, dubar, du); % Equation 6 
G	= compute_G(p, u, wp); % Equation 7 

% MEAN
rhs_ubar = ELu; 

% SPATIAL BASIS
rhs_u = (p - u*G)/C;

% STOCHASTIC COEFFICIENT
rhs_Y = h;

end
