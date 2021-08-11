%% Version 1: without integration by parts
function ADV = gen_advection_matrix(np, fe_mesh)
% This function returns the time-independent part of matrix system arising 
% from the spectral element discretization of the advection diffusion
% equation.
% 
% case (1): du/dt = a(x)*du/dx
%       MM*(u_new-u_old)/dt = a.*(ADV*u_old)
% 
% where MM is the mass matrix, a(x) is known function of x.
% 
%
% Input:
%          np -> order of the spectral elements
%     fe_mesh -> finite element mesh
%
% Output: ADV -> global advection matrix

[~, ~, w, ~, D] = gen_global_coordinate_system(np, fe_mesh);
nFE = length(fe_mesh)-1;

ADV = zeros(nFE*np+1);
for i=1:nFE
    wl = w((i-1)*(np+1)+1:i*(np+1));
    A = diag(wl)*D;
    ADV((i-1)*np+1:i*np+1,(i-1)*np+1:i*np+1) = ADV((i-1)*np+1:i*np+1,...
        (i-1)*np+1:i*np+1)+A;
end

end


%% Version 2 integrate by parts
% function ADV = gen_advection_matrix(np, fe_mesh)
% % This function returns the time-independent part of matrix system arising 
% % from the spectral element discretization of the advection diffusion
% % equation.
% % 
% % case (1): du/dt = a(x)*du/dx
% %       u_new = u_old + dt*MM\(ADV'*(a.*u_old)+ADV*a.*u_old)
% % 
% % case (2): du/dt = d(a(x)*u)/dx
% %       u_new = u_old + dt*MM\ADV'*(a.*u_old)
% % 
% % where MM is the mass matrix, a(x) is known function of x.
% % 
% %
% % Input:
% %          np -> order of the spectral elements
% %     fe_mesh -> finite element mesh
% %
% % Output: ADV -> global advection matrix
% 
% [~, ~, w, ~, D] = gen_global_coordinate_system(np, fe_mesh);
% nFE = length(fe_mesh)-1;
% 
% B = zeros(nFE*np+1);
% for i=1:nFE
%     wl = w((i-1)*(np+1)+1:i*(np+1));
%     A = diag(wl)*D;
%     B((i-1)*np+1:i*np+1,(i-1)*np+1:i*np+1) = B((i-1)*np+1:i*np+1,(i-1)*np+1:i*np+1)+A;
% end
% 
% % case (1)
% B(1,1) = B(1,1)+0.5;
% B(end,end) = B(end,end)-0.5;
% ADV = -B;
% 
% % case (2)
% % B(1,1) = B(1,1)+1;
% % B(end,end) = B(end,end)-1;
% % ADV = -B;
% 
% end