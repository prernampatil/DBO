% Code to solve the Stochastic Linear Advection equation
% Comparing results from 1) Analytical solution 2) DO  and 3) DBO
% Code modified for DBO by: Prerna Patil
% Original DO Code        : Hessam Babaee
% Written on              : 18th Jan 2019
% Last modified           : 12th September 2019

clc
clear
close all
clear global

%% Problem set-up
LW = 'linewidth';
RGB1 = [0 113 188]/norm([0 113 188]);
RGB2 = [216 82 24]/norm([216 82 24]);
RGB3 = [20 82 24 ]/norm([20 82 24 ]);

global xx
global du dubar dY x wp xr wr C Nr
global bY bu bubar bSigma
global vmean
global sigma

N    = 2;           % # of the basis {u_i}
Ns   = 128;         % the number of collocation points in physical space
Nr   = 64;          % the number of collocation points in random space
L    = 2*pi;        % the length of the interval [0,1]
rDim = 1;           % dimension of random space
xx   = chebfun('xx',[0,L]);
dom  = [0, L];

% Using uniform distibution for RVs between [-1,1]
vmean  = 1.0;
sigma  = 1.0;

x  = L*(0:Ns-1)'/Ns;   % collocation points
wp = L/Ns*ones(length(x),1);

tf = 10;            % final time
dt = 0.005;
t0 = 0+dt;
nTimeStep = ceil((tf-t0) / dt);

% Using uniform distribution between [-1,1]
[xr, wr] = lgwt(Nr, -1, 1);  wr = wr/2;

% Stochastic velocity
% c= vmean + sigma*xr;
% Evaluated at every legendre quadrature point xr(i)

% Initialization:
% @time t=0: The fluctuating parts are zero. Hence Y and u_i cannot be
% declared. Hence we start the simulation from 0+dt

% Initial condition for DO:
% Analytical solution: (Both computational and analytical solution use
% analytical solution for IC)
dubar  = u_exact(t0);
dubar  = dubar(x,:);
dY     = YDO_exact(t0);
du     = UDO_exact(t0);  du = du(x,:);
ndY    = dY;
ndu    = du;
ndubar = dubar;

% Initial conditions for DBO
% The solutions are obtained by using tranformations between DO->DBO
[nbu, nbSigma, nbY] = getDBO(du, dY, wr, wp);
bubar  = dubar;
bY     = nbY;
bu     = nbu;
bSigma = nbSigma;

% Declaring variables for error calculations
Error_meando  = zeros(1,nTimeStep);
Error_meandbo = zeros(1,nTimeStep);
Error_vardo   = zeros(1,nTimeStep);
Error_vardbo  = zeros(1,nTimeStep);
tic;
n=1;
t3=t0;
%Time loop:
while n <= nTimeStep

    % DO computation:
    % RK3 time stepping:
    t1 = dt * (n-1)+t0;
    [rhs_dubar1, rhs_du1, rhs_dY1] = compute_rhs_do(dubar, du, dY, N, Ns, Nr, xr, wr, wp, t1);
    
    t2 = dt * (n-0.5)+t0;
    dubar2 	= dubar + dt*rhs_dubar1/2.0;
    du2		= du    + dt*rhs_du1/2.0;
    dY2		= dY    + dt*rhs_dY1/2.0;
    [rhs_dubar2, rhs_du2 ,rhs_dY2] = compute_rhs_do(dubar2, du2, dY2, N, Ns, Nr, xr, wr, wp, t2);
    
    t3 = dt * n+t0;
    dubar3 	= dubar - dt*rhs_dubar1 + 2.0*dt*rhs_dubar2;
    du3		= du    - dt*rhs_du1    + 2.0*dt*rhs_du2;
    dY3		= dY    - dt*rhs_dY1    + 2.0*dt*rhs_dY2;
    [rhs_dubar3, rhs_du3, rhs_dY3] = compute_rhs_do(dubar3, du3, dY3, N, Ns, Nr, xr, wr, wp, t3);
    
    ndubar 	= dubar + dt*(rhs_dubar1+4.0*rhs_dubar2+rhs_dubar3)/6.0;
    ndu		= du    + dt*(rhs_du1   +4.0*rhs_du2   +rhs_du3)/6.0;
    ndY		= dY    + dt*(rhs_dY1   +4.0*rhs_dY2   +rhs_dY3)/6.0;
    
    C = ComputeCovBasis(ndY,wr);
    
    % DBO computations:
    % RK3 time stepping:
    t1 = dt * (n-1)+t0;
    [rhs_bubar1, rhs_bu1, rhs_bY1, rhs_bSigma1] = compute_rhs_dbo(bubar, bu, bY, bSigma, xr, wr, wp, t1);
    
    t2 = dt * (n-0.5)+t0;
    bubar2  = bubar  + dt*rhs_bubar1/2.0;
    bu2     = bu     + dt*rhs_bu1/2.0;
    bY2     = bY     + dt*rhs_bY1/2.0;
    bSigma2	= bSigma + dt*rhs_bSigma1/2.0;
    [rhs_bubar2, rhs_bu2, rhs_bY2, rhs_bSigma2] = compute_rhs_dbo(bubar2, bu2, bY2, bSigma2,xr, wr, wp, t2);
    
    t3 = dt * n+t0;
    bubar3 	= bubar  - dt*rhs_bubar1  + 2.0*dt*rhs_bubar2;
    bu3		= bu     - dt*rhs_bu1     + 2.0*dt*rhs_bu2;
    bY3		= bY     - dt*rhs_bY1     + 2.0*dt*rhs_bY2;
    bSigma3 = bSigma - dt*rhs_bSigma1 + 2.0*dt*rhs_bSigma2;
    [rhs_bubar3, rhs_bu3, rhs_bY3, rhs_bSigma3] = compute_rhs_dbo(bubar3 , bu3, bY3, bSigma3,xr, wr, wp, t3);
    
    nbubar 	= bubar  + dt*(rhs_bubar1 + 4.0*rhs_bubar2  +rhs_bubar3)/6.0;
    nbu		= bu     + dt*(rhs_bu1    + 4.0*rhs_bu2     +rhs_bu3)/6.0;
    nbY		= bY     + dt*(rhs_bY1    + 4.0*rhs_bY2     +rhs_bY3)/6.0;
    nbSigma = bSigma + dt*(rhs_bSigma1+ 4.0*rhs_bSigma2 +rhs_bSigma3)/6.0;
    
    % Enforce zero mean condition
    for i=1:N
        nbY(:,i) = nbY(:,i) - sum(nbY(:,i).*wr);
    end
    
    % Enforce gram schmidt condition on bY and bu
    [nbY, nbu] = gramSchmidt(nbY, nbu);
    
    cov_dbo = eig(nbSigma*nbSigma');
    
    % Exact solutions for DO equations:
    uedo  = u_exact(t3);
    uedo  = uedo(x);
    duedo = UDO_exact(t3);  duedo = duedo(x,:);
    dYedo = YDO_exact(t3);
    
    % Exact solutions for the DBO equations:
    % Obtained by using the transformations DO->DBO
    [duedbo, dSigmaedbo, dYedbo] = getDBO(duedo, dYedo, wr, wp);
    uedbo = uedo;
    
    % Calculate the error
    % L_2 norm for mean error
    err_meando  = ndubar - uedo;
    err_meandbo = nbubar - uedbo;
        
    Error_meando(n) = sqrt(sum(err_meando.*err_meando.*wp));
    Error_meandbo(n) = sqrt(sum(err_meandbo.*err_meandbo.*wp));
    
    % Calculate the weighted Frobenius norm for variance
    err_var1   = ndu*ndY' - duedo*dYedo';
    err_var2   = nbu*nbSigma*nbY' - duedbo*dSigmaedbo*dYedbo';
    err_vardo  = err_var1.*err_var1;
    err_vardbo = err_var2.*err_var2;
    
    Error_vardo(n)  = sqrt(wp'*err_vardo*wr);
    Error_vardbo(n) = sqrt(wp'*err_vardbo*wr);
    
    % DO update
    dY    = ndY;
    dubar = ndubar;
    du    = ndu;
    
    % DBO update
    bY     = nbY;
    bubar  = nbubar;
    bu     = nbu;
    bSigma = nbSigma;
    
    if mod(n,100)==0
        disp(['t=' num2str(n*dt) ' is being processed'])
    end
    
    n=n+1;	% necessary for while statement
    
end
toc
% Plot error for mean, variance and eigenvalues
% The figures are saved in Folder: Figures
if(~isfolder('ErrorPlots'))
    mkdir('ErrorPlots');
end
Time = t0:dt:tf-dt;

figure(1)
semilogy(Time(1:5:end),Error_meando(1:5:end),':','color',RGB2, LW, 2.5)
hold on
semilogy(Time(1:5:end),Error_meandbo(1:5:end),'.-','color',RGB1, LW,1.5)
hold on
title(['Error Mean t=',num2str(t3)])
xlim([t0 tf])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
xlabel('Time')
ylabel('L_2 error')
legend('DO', 'DBO')
drawnow
saveas(gcf,'ErrorPlots/MeanError','epsc')

figure(2)
semilogy(Time(1:5:end),Error_vardo(1:5:end),':','color',RGB2, LW, 2.5)
hold on
semilogy(Time(1:5:end),Error_vardbo(1:5:end),'.-','color',RGB1, LW, 1.5)
xlim([t0 tf])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
xlabel('Time')
ylabel('L_2 error')
legend('DO', 'DBO')
title(['Error variance t=',num2str(t3)])
drawnow
saveas(gcf,'ErrorPlots/VarError','epsc')

figure(3)
% Computing the analytical eigenvalues
tt = linspace(t0,tf,nTimeStep);
Lambda(:,2) = pi/2*(1 - sin(2*sigma*tt*pi)./(2*sigma*tt*pi));
Lambda(:,1) = pi/2*(1 + sin(2*sigma*tt*pi)./(2*sigma*tt*pi) - 2*(sin(tt*sigma*pi)).^2./(tt*sigma*pi).^2);
plot(tt,Lambda, LW, 1.5);
drawnow
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
title('Eigenvalues')
xlabel('Time')
ylabel('Eigenvalues')
drawnow
saveas(gcf,'ErrorPlots/Eigenvalues','epsc')