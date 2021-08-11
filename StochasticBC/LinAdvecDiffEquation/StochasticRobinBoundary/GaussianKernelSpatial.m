% Gaussian kernel for the Initial conditions 
% Written by: Prerna Patil 
% Date: 26th April 2021 

clc
clear 
clear global
close all 

%% -------------- Space Discretization ---------------
np                 = 4;
nFE                = 101;
Xmin               = 0; Xmax = 5;
fe_mesh            = linspace(Xmin, Xmax, nFE+1);
[~, x, we, wp, De] = gen_global_coordinate_system(np, fe_mesh);
we                 = we(1:np+1);
N                  = length(x);

% Form the kernel 
[x1, x2] = meshgrid(x,x);
lc       = 1; 
K        = (exp(-0.5*(x1-x2).^2/lc^2));
[Ux,Lx]   = eig(K*diag(wp)); Ux = real(Ux); Lx = real(Lx);
Lx        = diag(Lx);
[Lx, I]   = sort(Lx,'descend');
plot(cumsum(Lx)/sum(Lx),'*k','MarkerSize',1.5);
xlim([1,100]);
% Consider 99% energy of the modes 
Neigs = min(find(cumsum(Lx)/sum(Lx)>0.9999));
disp(Neigs)
Lx = Lx(1:Neigs);
Ux = Ux(:,I(1:Neigs));
Ux = Ux*diag(1./sqrt(wp'*(Ux.^2)));

tempstr = sprintf('SpatialGaussKern_d%d.mat',Neigs);
save(tempstr,'Lx','Ux');


