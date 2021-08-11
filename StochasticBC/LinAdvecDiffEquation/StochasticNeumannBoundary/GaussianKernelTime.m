% Gaussian kernel (for time)
% Written by: Prerna Patil 
clear
clear global 
close all 
clc

Tf = 5; 
dt = 0.0005; 
Ts = dt;
t  = 0:dt:Tf;
Nt  = size(t,2);
% Form the kernel 
wy = ones(Nt,1)*1/Nt;
[t1, t2] = meshgrid(t,t);
lc = 1.0; 
K  = (exp(-0.5*(t1-t2).^2/lc^2));
[Ut,L] = eig(K*diag(wy)); Ut = real(Ut); L = real(L);
L      = diag(L);
[L, I] = sort(L,'descend');
plot(cumsum(L)/sum(L),'*k','MarkerSize',1.5);
xlim([1,100]);
% Consider 99% energy of the modes 
Neigs = min(find(cumsum(L)/sum(L)>0.9999));
disp(Neigs)
L = L(1:Neigs);
Ut = Ut(:,I(1:Neigs));
Ut = Ut*diag(1./sqrt(wy'*(Ut.^2)));

tempstr = sprintf('GaussKern_d%d_Tf%d.mat',Neigs,Tf);
save(tempstr,'L','Ut');
