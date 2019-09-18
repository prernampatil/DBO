% Solution computed for PCM and stored for error analysis
% Modified by: Prerna Patil
% Original Code: Dr. Hessam Babaee
% Date: 6th Feb 2019
% Legnedre collocation points in random space
% u_t + u u_x = nu u_xx + (1+\xi)/2 * sin(2*pi*t)
% where
% u(x,0) = 0.5*(exp(cos(2*x))-1.5) *sin(3*x+2*pi*0.37);
% u(0,t) = u(2*pi,t)  % Periodic BC
% \xi \in U[-1, 1]
% 4th order RK is used as time integrator


clc
clear
close all

Ns = 128;   % Number of collocation points in physical space
Nr = 64;    % Numnber of points in random space
L = 2*pi;   % Length of the domain
x = L*(0:Ns-1)/Ns;
wp= L/Ns*ones(length(x),1);

nu =0.04;	% Diffusion coefficient
tf = 3;
t0 = 2;
dt = 0.001;
Nstep  = ceil(tf/dt);
NiniStep = ceil(t0/dt);
% For the KL solution
nEig = 9; % Number of modes computed for the KL solution

%Initial condition:
u0 = 0.5*(exp(cos(2*x))-1.5).*sin(3*x+2*pi*0.37);
[xr, wr] = lgwt(Nr, -1, 1);  wr = wr/sum(wr);

% Declare the matrices for the solution
u    = zeros(Ns, Nstep+1, Nr);      % Complete solution
mean = zeros(Ns, Nstep+1);          % Mean solution
var  = zeros(Ns, Nstep+1);          % Variance of the solution
tic
for i=1:Nr
    u(:,:,i) = burgers_solver_rk4(Ns, x, u0, tf, dt, nu, xr(i));
    mean = mean + u(:,:,i)*wr(i);
    var = var + u(:,:,i).^2* wr(i);
    if(mod(i,50) == 0)
        disp([num2str(i) ' is processed...']);
    end
    
    disp([num2str(i) ' completed!']);
    
end
toc
tic 
var = var - mean.^2;
% Remove data from the initial time time steps
nrem  = t0/dt;
mean  = mean(:,nrem+1:end);
var   = var(:,nrem+1:end);
u     = u(:,nrem+1:end,:);
Nstep = Nstep - nrem;
%Compute the KL solution at each time step 
for i=1:Nstep
    [upcm, ypcm] = getKL(reshape(u(:,i,:),Ns, Nr), reshape(mean(:,i),Ns,1), nEig, wr, wp);
    cov = eig(ComputeCovBasis(ypcm,wr));
    cov_pcm(:,i) = sort(cov,'descend');
end
fns = sprintf('pcm_rk4_Ns%d_Nr%d_nu%f_N%d_tf_%f.mat', Ns, Nr, nu, nEig, tf );% Store the solution
save(fns, 'mean', 'var', 'cov_pcm', 'u')
toc