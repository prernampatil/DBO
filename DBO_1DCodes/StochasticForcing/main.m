% Burger's equation random forcing
% Original code: Dr. Hessam Babaee
% Adapted for DBO by: Prerna Patil
% Date: 15th Jan 2019
% Code to evaluate stochastic forcing on the burgers equation
% The results are compared from three methods:
% 1)PCM    2)DO   3)DBO

%*******************************************************%
% Note: Run PCMrun/PCMmain.m with the appropriate       %
% parameters before running this code                   %
%*******************************************************%
clc
clear
close all
clear global
set(0,'defaultTextInterpreter','none');

% Initialization
%% Problem set-up
LW = 'linewidth';
RGB1 = [0 113 188 ]/norm([0 113 188 ]);
RGB2 = [216 82 24 ]/norm([216 82 24 ]);
RGB3 = [20 82 24  ]/norm([20  82 24 ]);

N   = 9;         % # of the basis {u_i}
Ns  = 128;       % the number of collocation points in physical space
Nr  = 64;        % the number of collocation points in random space
M   = Nr;
L   = 2*pi;      % the length of the interval [0,2*pi]
rDim= 1;         % dimension of random space
dom = [0, L];
mu  = 0.04;      % diffusion coefficient

x   = L*(0:Ns-1)'/Ns;   % collocation points
wp  = L/Ns*ones(length(x),1);
t0  = 0;
tf  = 3;                % final time
dt  = 0.001;

[xr, wr] = lgwt(Nr, -1, 1);  wr = wr/sum(wr);
%Initial condition:
u0 = 0.5*(exp(cos(2*x))-1.5).*sin(3*x+2*pi*0.37);
ts = 2;
nTimeStep  = ceil((tf-ts)/dt);

% Calculate the solution using PCM till the stochasticity sets in:
disp('Computing PCM solution till switching time')
[u_pcm, mean, var]= initialPCM(Ns, Nr, xr, x, u0,t0, dt, mu, wr,ts); %Changing the ICs for the cases
disp('Solution computed till switching time')
u0 = mean;
t0 = ts;
u_pcm0 = u_pcm;
tstart = t0;
% Get KL:
[ndu, ndY] = getKL(u_pcm, mean, N, wr, wp);
ndubar = mean;
dubar  = mean;
dY  = ndY;
du  = ndu;
PCM0 =  eig(ComputeCovBasis(ndY,wr));
% Get DBO decompostion:
[nbu, nbSigma, nbY] = getDBO(ndu, ndY, wr, wp);
nbubar = mean;
bubar  = mean;
bY = nbY;
bu = nbu;
bSigma = nbSigma;

n   = 1;
t3  = t0;
tic;
DO  = true;
DBO = true;
PCM = true;
fns = sprintf('PCMRun/pcm_rk4_Ns%d_Nr%d_nu%f_N%d_tf_%f.mat', Ns, Nr, mu, N, tf);
load(fns);
while(n <= (nTimeStep-1))
    
    % Compute DO solutions:
    if(DO)
        cov_dy = eig(ComputeCovBasis(dY,wr));
        
        t1 = dt * (n-1)  + t0;
        [rhs_dubar1, rhs_du1, rhs_dY1, du,  dY]  = compute_rhs_do(dubar, du, dY, N, Ns, Nr, xr, wr, wp, t1, mu);
        
        t2 = dt * (n-1/2)+ t0;
        dubar2 	= dubar  + dt*rhs_dubar1/2.0;
        du2		= du     + dt*rhs_du1/2.0;
        dY2		= dY     + dt*rhs_dY1/2.0;
        [rhs_dubar2, rhs_du2, rhs_dY2, du2, dY2] = compute_rhs_do(dubar2, du2, dY2, N, Ns, Nr, xr, wr, wp, t2, mu);
        
        t3 = dt * n+t0;
        dubar3 	= dubar  - dt*rhs_dubar1 + 2.0*dt*rhs_dubar2;
        du3		= du     - dt*rhs_du1    + 2.0*dt*rhs_du2;
        dY3		= dY     - dt*rhs_dY1    + 2.0*dt*rhs_dY2;
        [rhs_dubar3, rhs_du3, rhs_dY3, du3, dY3] = compute_rhs_do(dubar3, du3, dY3, N, Ns, Nr, xr, wr, wp, t3, mu);
        
        ndubar 	= dubar  + dt*(rhs_dubar1+ 4.0*rhs_dubar2+rhs_dubar3)/6.0;
        ndu		= du     + dt*(rhs_du1   + 4.0*rhs_du2   +rhs_du3)/6.0;
        ndY		= dY     + dt*(rhs_dY1   + 4.0*rhs_dY2   +rhs_dY3)/6.0;
        
        C=ComputeCovBasis(ndY,wr);
    end
    
    % Compute DBO solutions:
    if(DBO)
        t1 = dt * (n-1)+t0;
        [rhs_bubar1, rhs_bu1, rhs_bY1, rhs_bSigma1] = compute_rhs_dbo(bubar, bu, bY, bSigma, mu,xr, wr, wp, t1);
        
        t2 = dt * (n-1/2)+t0;
        bubar2  = bubar  + dt*rhs_bubar1/2.0;
        bu2     = bu     + dt*rhs_bu1/2.0;
        bY2     = bY     + dt*rhs_bY1/2.0;
        bSigma2	= bSigma + dt*rhs_bSigma1/2.0;
        [rhs_bubar2, rhs_bu2, rhs_bY2, rhs_bSigma2] = compute_rhs_dbo(bubar2, bu2, bY2, bSigma2, mu,xr, wr, wp, t2);
        
        t3 = dt * n+t0;
        bubar3 	= bubar  - dt*rhs_bubar1  + 2.0*dt*rhs_bubar2;
        bu3		= bu     - dt*rhs_bu1     + 2.0*dt*rhs_bu2;
        bY3		= bY     - dt*rhs_bY1     + 2.0*dt*rhs_bY2;
        bSigma3 = bSigma - dt*rhs_bSigma1 + 2.0*dt*rhs_bSigma2;
        [rhs_bubar3, rhs_bu3, rhs_bY3, rhs_bSigma3] = compute_rhs_dbo(bubar3 , bu3, bY3, bSigma3, mu,xr, wr, wp, t3);
        
        nbubar 	= bubar  + dt*(rhs_bubar1 + 4.0*rhs_bubar2  +rhs_bubar3)/6.0;
        nbu		= bu     + dt*(rhs_bu1    + 4.0*rhs_bu2     +rhs_bu3)/6.0;
        nbY		= bY     + dt*(rhs_bY1    + 4.0*rhs_bY2     +rhs_bY3)/6.0;
        nbSigma = bSigma + dt*(rhs_bSigma1+ 4.0*rhs_bSigma2 +rhs_bSigma3)/6.0;
        % Enforce zero mean condition
        for i=1:N
            nbY(:,i) = nbY(:,i) - sum(nbY(:,i).*wr);
        end
        
        % Enforce gram schmidt condition on bY and bu
        nbY(:,1) = nbY(:,1)/sum(nbY(:,1).*nbY(:,1).*wr);
        for i=2:N
            tempY = nbY(:,i);
            tempU = nbu(:,i);
            for j=1:i-1
                tempY = tempY - sum( nbY(:,i).*nbY(:,j).*wr )/sum( nbY(:,j).*nbY(:,j).*wr )*nbY(:,j);
                tempU = tempU - sum( nbu(:,i).*nbu(:,j).*wp )/sum( nbu(:,j).*nbu(:,j).*wp )*nbu(:,j);
            end
            tempY = tempY/ sum(tempY.*tempY.*wr);
            nbY(:,i) = tempY;
            tempU = tempU/ sum(tempU.*tempU.*wp);
            nbu(:,i) = tempU;
        end
        cov_dbo = eig(nbSigma*nbSigma');
    end
    
    % Compute error:
    if(DO)
        err_meando  = ndubar - mean(:,n+1);
        err_var1    = ndu*ndY'- reshape(u(:,n+1,:),Ns,Nr)+ mean(:,n+1);
        err_vardo   = err_var1.*err_var1;
        ErrorDO(n)     = sqrt(err_meando'*(wp.*err_meando));
        VarErrDO(n)    = sqrt(wp'*err_vardo*wr);
        LambdaDO(:,n) = cov_dy;
    end
    if(DBO)
        err_meandbo = nbubar - mean(:,n+1);
        err_var2    = nbu*nbSigma*nbY'- reshape(u(:,n+1,:),Ns,Nr)+ mean(:,n+1);
        err_vardbo  = err_var2.*err_var2;
        ErrorDBO(n)    = sqrt(err_meandbo'*(wp.*err_meandbo));
        VarErrDBO(n)   = sqrt(wp'*err_vardbo*wr);
        LambdaDBO(:,n) = cov_dbo;
    end
    
    % DO update
    if(DO)
        dY = ndY;
        dubar = ndubar;
        du = ndu;
    end
    % DBO update
    if(DBO)
        bY     = nbY;
        bubar  = nbubar;
        bu     = nbu;
        bSigma = nbSigma;
    end
    if mod(n,100)==0
        drawnow
        disp(['t=' num2str(ts+n*dt) ' is being processed'])
    end
    n=n+1; % necessary for while statement
    
end
% Plot error for mean, variance and eigenvalues
% The figures are saved in Folder: Figures
if(~isfolder('ErrorPlots'))
    mkdir('ErrorPlots');
end
Time = t0:dt:tf-dt-dt;
figure(1)
if(DBO)
    semilogy(Time,VarErrDBO,'--','color', RGB1, LW, 1.5);
    hold on
end
if(DO)
    semilogy(Time,VarErrDO,':','color',RGB2, LW,2.5 );
    hold on
end
xlabel('Time')
ylabel('$\mathrm{L_2}$ ','interpreter', 'latex')
legend('DBO','DO')
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
saveas(gcf,'ErrorPlots/VarError','epsc')

figure(2)
if(DBO)
    semilogy(Time,ErrorDBO,'--','color', RGB1, LW, 1.5);
    hold on
end
if(DO)
    semilogy(Time,ErrorDO,':','color',RGB2, LW, 2.5 );
    hold on
end
xlabel('Time')
ylabel('$\mathrm{L_2}$ Error','Interpreter','latex')
legend('DBO','DO')
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
saveas(gcf,'ErrorPlots/MeanError','epsc')

figure(3)
if(PCM)
    semilogy(Time, cov_pcm(:,1:end-1), '-k', LW, 1.5);
    hold on
end
if(DBO)
    semilogy(Time, LambdaDBO, '--', 'color', RGB1, LW, 1.5);
    hold on
end
if(DO)
    semilogy(Time, LambdaDO,':', 'color', RGB2, LW, 2.5);
    hold on
end
h = zeros(3,1);
h(1) = plot(NaN, NaN, ':', 'color', RGB2, LW, 2.5);
h(2) = plot(NaN, NaN, '.-', 'color', RGB1, LW, 1.5);
h(3) = plot(NaN, NaN, '-k', LW, 1.5);
legend(h,'DO','DBO','PCM');
xlabel('Time')
ylabel('Eigenvalues')
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
saveas(gcf,'ErrorPlots/Eigenvalues','epsc')