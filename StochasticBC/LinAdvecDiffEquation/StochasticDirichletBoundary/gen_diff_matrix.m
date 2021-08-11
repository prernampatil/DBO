% % This function returns the time-independent part of matrix system arising
% % from the spectral element discretization of the diffusion equation L(u) = f(x)
% % L(u) = -d(a(x)*du/dx)/dx
% % where a(x) is known function of x.
% %
% % Weak formulation:
% % (L(u), v) = (f(x), v) ==> DIF*u = w.*f(x)
% % where (.,.) means numerical quadrature rule
% %
% % Input:
% %          np -> order of the spectral elements
% %         nFE -> No. of finite elements
% %           w -> quadrature weights of integration over one single element
% %           D -> element differentiation matrix
% %           a -> vector of diffusion coefficients
% %
% % Output: DIF -> global advection matrix
%% Version 1 (Correct)
% First write down the bi-linear form of weak formulation, then do projection
function DIF = gen_diff_matrix(np, nFE, w, D, a, bc_left, bc_right)
DIF = zeros(nFE*np+1);
for i=1:nFE
    %     al = a((i-1)*np+1:i*np+1);
    %     A = D'*diag(w.*al)*D;
    A = D'*diag(w*a)*D;
    DIF((i-1)*np+1:i*np+1,(i-1)*np+1:i*np+1) = ...
        DIF((i-1)*np+1:i*np+1,(i-1)*np+1:i*np+1)+A;
end

% if periodic boundary condition, no need for these two lines
if (bc_left =='D')
    DIF(1, 1:np+1) = DIF(1, 1:np+1)+a(1)*D(1,:);
end
if (bc_right =='D')
    DIF(end, end-np:end) = DIF(end, end-np:end)-a(end)*D(end,:);
end



end

%% Version 2 (Wrong)
% First do projection, then write as a bi-linear form of weak formulation
% function DIF = gen_diff_matrix(np, nFE, w, D, a)
%
% DIF = zeros(nFE*np+1);
% for i = 1:nFE
%     index = (i-1)*np+1:i*np+1;
%     al = a(index);
%     A = D'*diag(w.*al)*D;
%     DIF(index,index) = DIF(index,index) + A;
%     DIF(index(1),   index) = DIF(index(1),   index) + al(1)*D(1,:);
%     DIF(index(end), index) = DIF(index(end), index) - al(end)*D(end,:);
% end
% end