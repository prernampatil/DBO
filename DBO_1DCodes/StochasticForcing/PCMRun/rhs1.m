% u_n(x,t) = \sum_{n=0}^{N-1} \phi_n(t) h_n(x)
%
% \phi = (\phi_0(t), \phi_1(t), ..., \phi_{N-1}(t)
% 
% Inserting u_n into the Burger's equation, we get system of ODE for
% \phi.
%
%   d\phi/dt = - diag(\phi) * D*\phi + \nu * D^2*\phi
% where 
%   D is the Fourier differentiation matrix.
%
% We compute D^m*\phi through FFT.
%

function r = rhs1(u, N, nu, xi, t)

%    d1 = FourierDerivativeByFFT(N, u, 1);
%    d2 = FourierDerivativeByFFT(N, u, 2);        

   d1 = fourdifft(u, 1);
   d2 = fourdifft(u, 2);
   r = -u.*d1 + nu*d2 + 0.5*(xi+1)*sin(2*pi*t);     
end
