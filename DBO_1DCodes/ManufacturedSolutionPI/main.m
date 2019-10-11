% Code to solve the Stochastic Burgers' equation
% Results are compared for 1.) DO 2.) PI-DO 3.) DBO 4.) Analytical solution
% Code modified for DBO by: Prerna Patil
% Original DO Code        : Hessam Babaee
% Written on              : 18th Jan 2019
% Last modified           : 13th September 2019

clc
clear
close all
clear global

%% Problem set-up
LW = 'linewidth';
global RGB1 RGB2 RGB3 RGB4 RGB5
RGB1 = [0  113 188]/norm([0 113 188]);
RGB2 = [216 82 24 ]/norm([216 82 24]);
RGB3 = [0 153 153 ]/norm([0 153 153]);
RGB4 = [38, 115, 77]/norm([38, 115, 77]);
RGB5 = [122, 0, 153]/norm([122, 0, 153]);

global w_n s1x c1x s2x c2x
global Y du dupi pu dubar dubarpi pubar x wp xr wr nu C Nr
global epsilon epname
global tol_PI
tol_PI = 10^-10;

epsilon = 10^-3;
if(epsilon == 10^-5)
    epname = 'epsilon10_5';
elseif(epsilon == 10^-3)
    epname = 'epsilon10_3';
end
videoOutputrequired = 1;
if(videoOutputrequired)
    % Create video
    writerObj1 = VideoWriter('Phasespace.avi');
    writerObj1.FrameRate = 30;
    writerObj1.Quality = 80;
    open(writerObj1);
    writerObj2 = VideoWriter('SpatialModes.avi');
    writerObj2.FrameRate = 30;
    writerObj2.Quality = 80;
    open(writerObj2);
end
Snapshots = [200, 1200, 3200, 5200];
% Create folder for saving the results
CreateFolders()

xx     = chebfun('xx',[0,pi]);
N      = 2;             % # of the basis {u_i}
Ns     = 128;           % the number of collocation points in physical space
Nr     = 1024;          % the number of collocation points in random space
L      = 2*pi;          % the length of the interval [0,2*pi]
rDim   = 2;             % dimension of random space

dom =[0, L];

w_n=1; s1x =chebfun(@SX, dom,'trig');
w_n=2; s2x =chebfun(@SX, dom,'trig');

w_n=1; c1x =chebfun(@CX, dom,'trig');
w_n=2; c2x =chebfun(@CX, dom,'trig');

nu = 0.05;               % diffusion coefficient

x   = L*(0:Ns-1)'/Ns;    % collocation points
wp  = L/Ns*ones(length(x),1);
t0  = 0.01;
tf  = 3*pi;              % final time
dt  = 0.001;
nTimeStep  = ceil((tf-t0) / dt);

% Loading the stochastic collocation points and weights
% Calculated from the Me-PCM code
xr =  load('ColPnts_d2-1024.dat');
wr =  load('ColWgts_d2-1024.dat');

% Initialization for DO
dubar = u_exact(t0);       dubar = dubar(x,:);
dY    = YDO_exact(t0);
du    = UDO_exact(t0);     du = du(x,:);

ndY    = zeros(Nr, N);
ndu    = zeros(Ns, N);
ndubar = zeros(Ns, 1);

% Initialization for PI-DO
dubarpi = u_exact(t0);       dubarpi = dubarpi(x,:);
dYpi    = YDO_exact(t0);
dupi    = UDO_exact(t0);     dupi = dupi(x,:);

ndYpi    = zeros(Nr, N);
ndupi    = zeros(Ns, N);
ndubarpi = zeros(Ns, 1);

% Intialize for DBO
[npu, npSigma, npY] = getDBO(du, dY, wr, wp);
npubar = dubar;
pubar  = dubar;
pY     = npY;
pu     = npu;
pSigma = npSigma;

% Initialize for BO
[nbu, nbY] = getBO(du, dY, wr, wp);
nbubar = dubar;
bubar  = dubar;
bY     = nbY;
bu     = nbu;

% Initializing P and Q for the matrix differential eqns.
cov_bu = ComputeCovBasis(bu,wp);
inv_cov_bu = inv(cov_bu);
P = eye(N);
Q = diag(sqrt(diag(inv_cov_bu)));

tic;
n=1;

BO   = true;
DO   = true;
PIDO = true;
DBO  = true;
t3=0;

% Declare the variables for error
% DO
Errord = zeros(1, nTimeStep);
Error_vard = zeros(1, nTimeStep);
% PI-DO
Errordpi = zeros(1, nTimeStep);
Error_vardpi = zeros(1, nTimeStep);
%DBO
Error = zeros(1, nTimeStep);
Error_var = zeros(1, nTimeStep);
%BO
Errorp = zeros(1, nTimeStep);
Error_varp = zeros(1, nTimeStep);
% Declare variables for Eigenvalues
EigDO   = zeros(N, nTimeStep);
EigDOpi = zeros(N, nTimeStep);
EigBO   = zeros(N, nTimeStep);
EigDBO  = zeros(N, nTimeStep);
while n <= nTimeStep
    
    if(BO)
        cov_bu = ComputeCovBasis(bu,wp);
    end
    % Evolution equations for DO
    % RK3 time stepping
    if(DO)
        t1 = dt * (n-1)+t0;
        [rhs_dubar1, rhs_du1, rhs_dY1, du, dY]   = compute_rhs_do(dubar, du, dY, N, Ns, Nr, xr, wr, wp, t1, nu);
        
        t2 = dt * (n-1/2)+t0;
        dubar2 	= dubar + dt*rhs_dubar1/2.0;
        du2		= du    + dt*rhs_du1/2.0;
        dY2		= dY    + dt*rhs_dY1/2.0;
        [rhs_dubar2, rhs_du2, rhs_dY2, du2, dY2] = compute_rhs_do(dubar2, du2, dY2, N, Ns, Nr, xr, wr, wp, t2, nu);
        
        t3 = dt * n+t0;
        dubar3 	= dubar - dt*rhs_dubar1 + 2*dt*rhs_dubar2;
        du3		= du    - dt*rhs_du1    + 2*dt*rhs_du2;
        dY3		= dY    - dt*rhs_dY1    + 2*dt*rhs_dY2;
        [rhs_dubar3, rhs_du3, rhs_dY3, du3, dY3] = compute_rhs_do(dubar3, du3, dY3, N, Ns, Nr, xr, wr, wp, t3, nu);
        
        ndubar 	= dubar + dt*(rhs_dubar1+4*rhs_dubar2+rhs_dubar3)/6.0;
        ndu		= du    + dt*(rhs_du1   +4*rhs_du2   +rhs_du3)/6.0;
        ndY		= dY    + dt*(rhs_dY1   +4*rhs_dY2   +rhs_dY3)/6.0;
        
        C=ComputeCovBasis(ndY,wr);
        cov_dy = eig(C);
        clear rhs_dubar1 rhs_dubar2 rhs_dubar3
        clear rhs_du1 rhs_du2 rhs_du3
        clear rhs_dY1 rhs_dy2 rhs_dy3
        clear du2 du3 dY2 dY3 dubar2 dubar3
        % Enforcing the zero mean condition on the Y coeffs
        for i=1:N
            ndY(:,i) = ndY(:,i) - sum(ndY(:,i).*wr);
        end
    end
    % Evolution equations for PI-DO
    % RK3 time stepping
    if(PIDO)
        t1 = dt * (n-1)+t0;
        [rhs_dubarpi1, rhs_dupi1, rhs_dYpi1, dupi, dYpi]   = compute_rhs_dopi(dubarpi, dupi, dYpi, N, Ns, Nr, xr, wr, wp, t1, nu);
        
        t2 = dt * (n-1/2)+t0;
        dubarpi2 	= dubarpi + dt*rhs_dubarpi1/2.0;
        dupi2		= dupi    + dt*rhs_dupi1/2.0;
        dYpi2		= dYpi    + dt*rhs_dYpi1/2.0;
        [rhs_dubarpi2, rhs_dupi2, rhs_dYpi2, dupi2, dYpi2] = compute_rhs_dopi(dubarpi2, dupi2, dYpi2, N, Ns, Nr, xr, wr, wp, t2, nu);
        
        t3 = dt * n+t0;
        dubarpi3 	= dubarpi - dt*rhs_dubarpi1 + 2*dt*rhs_dubarpi2;
        dupi3		= dupi    - dt*rhs_dupi1    + 2*dt*rhs_dupi2;
        dYpi3		= dYpi    - dt*rhs_dYpi1    + 2*dt*rhs_dYpi2;
        [rhs_dubarpi3, rhs_dupi3, rhs_dYpi3, dupi3, dYpi3] = compute_rhs_dopi(dubarpi3, dupi3, dYpi3, N, Ns, Nr, xr, wr, wp, t3, nu);
        
        ndubarpi 	= dubarpi + dt*(rhs_dubarpi1+ 4*rhs_dubarpi2+ rhs_dubarpi3)/6.0;
        ndupi		= dupi    + dt*(rhs_dupi1   + 4*rhs_dupi2   + rhs_dupi3)/6.0;
        ndYpi		= dYpi    + dt*(rhs_dYpi1   + 4*rhs_dYpi2   + rhs_dYpi3)/6.0;
        
        Cpi=ComputeCovBasis(ndYpi,wr);
        [Cpi, EYYpi_inv, ndYpi] = pseudoinv(Cpi, ndYpi);
        cov_dypi = sort(eig(Cpi));
        clear rhs_dubarpi1 rhs_dubarpi2 rhs_dubarpi3
        clear rhs_dupi1 rhs_dupi2 rhs_dupi3
        clear rhs_dYpi1 rhs_dypi2 rhs_dypi3
        clear dupi2 dupi3 dYpi2 dYpi3 dubarpi2 dubarpi3
        % Enforcing the zero mean condition on the Y coeffs
        for i=1:N
            ndYpi(:,i) = ndYpi(:,i) - sum(ndYpi(:,i).*wr);
        end
    end
    % Evolution equations for BO
    % RK3 time stepping
    if(BO)
        t1 = dt * (n-1)+t0;
        [rhs_bubar1, rhs_bu1, rhs_bY1, bS] = compute_rhs_bo_v1(bubar, bu, bY, N, Ns, Nr, xr, wr, wp, t1, nu);
        diagS = diag(bS);
        
        t2 = dt * (n-1/2)+t0;
        bubar2 	= bubar + dt*rhs_bubar1/2.0;
        bu2		= bu    + dt*rhs_bu1/2.0;
        bY2		= bY    + dt*rhs_bY1/2.0;
        
        [rhs_bubar2, rhs_bu2, rhs_bY2, ~] = compute_rhs_bo_v1(bubar2, bu2, bY2, N, Ns, Nr, xr, wr, wp, t2, nu);
        
        t3 = dt * n+t0;
        bubar3 	= bubar - dt*rhs_bubar1 + 2*dt*rhs_bubar2;
        bu3		= bu    - dt*rhs_bu1    + 2*dt*rhs_bu2;
        bY3		= bY    - dt*rhs_bY1    + 2*dt*rhs_bY2;
        
        [rhs_bubar3, rhs_bu3, rhs_bY3, ~] = compute_rhs_bo_v1(bubar3, bu3, bY3, N, Ns, Nr, xr, wr, wp, t3, nu);
        
        nbubar 	= bubar + dt*(rhs_bubar1+4*rhs_bubar2 +rhs_bubar3)/6.0;
        nbu		= bu    + dt*(rhs_bu1   +4*rhs_bu2    +rhs_bu3)/6.0;
        nbY		= bY    + dt*(rhs_bY1   +4*rhs_bY2    +rhs_bY3)/6.0;
        clear rhs_bubar1 rhs_bubar2 rhs_bubar3
        clear rhs_bu1 rhs_bu2 rhs_bu3
        clear rhs_bY1 rhs_by2 rhs_by3
        clear bu2 bu3 bY2 bY3 bubar2 bubar3
        % Enforcing the zero mean condition on the Y coeffs
        for i=1:N
            nbY(:,i) = nbY(:,i) - sum(nbY(:,i).*wr);
        end
    end
    % Evolution equations for DBO
    % RK3 time stepping
    if(DBO)
        t1 = dt * (n-1)+t0;
        [rhs_pubar1, rhs_pu1, rhs_pY1, rhs_pSigma1] = compute_rhs_dbo(pubar, pu, pY, pSigma, nu,xr, wr, wp, t1);
        
        t2 = dt * (n-1/2)+t0;
        pubar2  = pubar  + dt*rhs_pubar1/2.0;
        pu2     = pu     + dt*rhs_pu1/2.0;
        pY2     = pY     + dt*rhs_pY1/2.0;
        pSigma2	= pSigma + dt*rhs_pSigma1/2.0;
        [rhs_pubar2, rhs_pu2, rhs_pY2, rhs_pSigma2] = compute_rhs_dbo(pubar2, pu2, pY2, pSigma2, nu,xr, wr, wp, t2);
        
        t3 = dt * n+t0;
        pubar3 	= pubar  - dt*rhs_pubar1  + 2.0*dt*rhs_pubar2;
        pu3		= pu     - dt*rhs_pu1     + 2.0*dt*rhs_pu2;
        pY3		= pY     - dt*rhs_pY1     + 2.0*dt*rhs_pY2;
        pSigma3 = pSigma - dt*rhs_pSigma1 + 2.0*dt*rhs_pSigma2;
        [rhs_pubar3, rhs_pu3, rhs_pY3, rhs_pSigma3] = compute_rhs_dbo(pubar3 , pu3, pY3, pSigma3, nu,xr, wr, wp, t3);
        
        npubar 	= pubar  + dt*(rhs_pubar1 + 4.0*rhs_pubar2  +rhs_pubar3)/6.0;
        npu		= pu     + dt*(rhs_pu1    + 4.0*rhs_pu2     +rhs_pu3)/6.0;
        npY		= pY     + dt*(rhs_pY1    + 4.0*rhs_pY2     +rhs_pY3)/6.0;
        npSigma = pSigma + dt*(rhs_pSigma1+ 4.0*rhs_pSigma2 +rhs_pSigma3)/6.0;
        clear rhs_pubar1 rhs_pubar2 rhs_pubar3
        clear rhs_pu1 rhs_pu2 rhs_pu3
        clear rhs_pY1 rhs_py2 rhs_py3
        clear pu2 pu3 pY2 pY3 pubar2 pubar3
        
        % Enforce zero mean condition
        for i=1:N
            npY(:,i) = npY(:,i) - sum(npY(:,i).*wr);
        end
        
        % Enforce gram schmidt condition on pY and pu
        [npY, npu] = gramSchmidt(npY, npu);
        cov_p = eig(npSigma*npSigma');
    end
    
    % Exact solution for DO
    % Using the analytical equations
    if(DO || PIDO)
        t3 = double(t3);
        ue = u_exact(t3);     ue = ue(x);
        due = UDO_exact(t3);  due = due(x,:);
        dYe = YDO_exact(t3);
    end
    % Exact solutions for the DBO equations:
    % Obtained using transformations DO->DBO
    if(DBO)
    [pue, pSigmae, pYe] = getDBO(due, dYe, wr, wp);
    end
    % Exact solutions for the BO equations:
    % Obtained using transformations DO->BO
    if(BO)
    [bue, bYe] = getBO(due, dYe, wr, wp);
    end
    % Calculate the error
    % L_2 norm for the mean
    % Frobenius norm for the variance
    if(DO)
        derr_mean     = ndubar - ue;
        Errord(n)     = sqrt(derr_mean'*(wp.*derr_mean));
        err_var1      = ndu*ndY' - due*dYe';
        derr_var      = err_var1.*err_var1;
        Error_vard(n) = sqrt(wp'*derr_var*wr);
        EigDO(:,n) = cov_dy;
    
    end
    if(PIDO)
        dpierr_mean     = ndubarpi - ue;
        Errordpi(n)     = sqrt(dpierr_mean'*(wp.*dpierr_mean));
        err_varpi1      = ndupi*ndYpi' - due*dYe';
        dpierr_var      = err_varpi1.*err_varpi1;
        Error_vardpi(n) = sqrt(wp'*dpierr_var*wr);
        EigDOpi(:,n)    = cov_dypi;
    
    end
    if(BO)
        err_mean      = nbubar - ue;
        Error(n)      = sqrt(err_mean'*(wp.*err_mean));
        err_var2      = nbu*nbY' - bue*bYe';
        err_var       = err_var2.*err_var2;
        Error_var(n)  = sqrt(wp'*err_var*wr);
        EigBO(:,n) = diag(cov_bu);
    
    end
    if(DBO)
        perr_mean     = npubar - ue;
        Errorp(n)     = sqrt(perr_mean'*(wp.*perr_mean));
        err_var3      = npu*npSigma*npY' - pue*pSigmae*pYe';
        perr_var      = err_var3.*err_var3;
        Error_varp(n) = sqrt(wp'*perr_var*wr);
        EigDBO(:,n) = cov_p;
    end
    
    % plot the spatial and stochastic basis
    if(mod(n,20)==0)
        % Compute the covaraince matrix for BO
        if(BO)
        C0 = zeros(N,N);
        for i=1:N
            for j=1:N
                C0(i,j) = sum(nbu(:,i) .* nbu(:,j) .* wp);
            end
        end
        [E,D,Et] = svd(C0);
        end
        figure(1)
        clf
        subplot(1,2,1)
        PlotSpatialModes(DO,PIDO, BO, DBO, D, E, nbu, ndu,ndupi, npu, pue,1);
        title('Physical space basis (Mode 1)')
        if(ismember(n,Snapshots))
            fname = sprintf('Basisplots/%s/PhyBasisMode1_%d',epname,n);
            print(fname,'-depsc');
        end
        subplot(1,2,2)
        PlotSpatialModes(DO,PIDO, BO, DBO, D, E, nbu, ndu,ndupi, npu, pue,2);
        title('Physical space basis (Mode 2)')
        set(gcf,'Position',[100 100 1000 500])
        if(ismember(n,Snapshots))
            fname = sprintf('Basisplots/%s/PhyBasisMode2_%d',epname,n);
            print(fname,'-depsc');
        end
        if(videoOutputrequired)
            frame = getframe(gcf);
            writeVideo(writerObj2,frame);
        end
        
        % Compute the reduced covariance matrix from DO
        if(DO)
        C0 = zeros(N,N);
        for i=1:N
            for j=1:N
                C0(i,j) = sum(ndY(:,i) .* ndY(:,j) .* wr);
            end
        end
        [E,D,Et] = svd(C0);
        YBO_1 = (inv(sqrt(D))*E'*ndY')'; % Converted from DO
        end
        if(DBO)
        YBO_2 = npY;
        end
        % Compute the reduced covariance matrix from DO
        if(PIDO)
        C0 = zeros(N,N);
        clear E D Et
        for i=1:N
            for j=1:N
                C0(i,j) = sum(ndYpi(:,i) .* ndYpi(:,j) .* wr);
            end
        end
        [E,D,Et] = svd(C0);
        YBOpi_1 = (inv(sqrt(D))*E'*ndYpi')'; % Converted from DO
        end
        
        figure(3)
        PlotPhasespace(DO, PIDO, BO, DBO,  YBO_1, YBOpi_1, nbY, YBO_2, bYe);
        if(ismember(n,Snapshots))
            fname = sprintf('Phasespace/%s/Phasespace_%d',epname,n);
            print(fname,'-depsc');
            saveas(gcf, fname, 'fig')
        end
        if(videoOutputrequired)
            frame = getframe(gcf);
            writeVideo(writerObj1,frame);
        end
    end
    
    clear ue bue bYe due dYe pue pSigmae pYe
    clear err_mean err_var err_meand err_vard perr_var perr_mean
    
    
    % update state
    if(BO)
        % BO update
        Y     = nbY;
        U     = nbu;
        bY    = nbY;
        bubar = nbubar;
        bu    = nbu;
        dubar = ndubar;
    end
    if(DO)
        % DO update
        Y     = ndY;
        U     = ndu;
        dY    = ndY;
        dubar = ndubar;
        du    = ndu;
    end
    if(PIDO)
        % DO update
        Ypi     = ndYpi;
        Upi     = ndupi;
        dYpi    = ndYpi;
        dubarpi = ndubarpi;
        dupi    = ndupi;
    end
    if(DBO)
        % DBO update
        pY     = npY;
        pubar  = npubar;
        pu     = npu;
        pSigma = npSigma;
    end
    
    if mod(n,100)==0
        disp(['t=' num2str(n*dt) ' is being processed'])
    end
    
    n=n+1;	% necessary for while statement
end
if(videoOutputrequired)
    % Clode the video files
    close(writerObj1);
    close(writerObj2);
    fname = sprintf('Phasespace/%s',epname);
    movefile('Phasespace.avi', fname);
    fname = sprintf('Basisplots/%s',epname);
    movefile('SpatialModes.avi', fname);
end
toc
figure(4)
% Calculate the eigenvalues
tt = linspace(t0,tf,nTimeStep);
Lambda(:,1) = (4.5 + sin(tt')).^2;
Lambda(:,2) = (epsilon*(1.5 + cos(3*tt'))).^2;
semilogy(tt, Lambda,'-k', LW, 1.5);

if(BO)
    semilogy(tt, EigBO, ':','color', RGB1, LW, 2.5);
    hold on
end
if(DO)
    semilogy(tt(1:400:end), EigDO(:,1:400:end),'d','color', RGB3,'MarkerSize',3.5);
    hold on
    semilogy(tt, EigDO,'-','color', RGB3, LW, 1.5);
end
if(PIDO)
    semilogy(tt(400:400:end), EigDOpi(:,400:400:end),'d','color', RGB5,'MarkerSize',3.5);
    hold on
    semilogy(tt, EigDOpi,'-.','color', RGB5, LW, 1.5);
end
if(DBO)
    semilogy(tt, EigDBO, '--','color', RGB2, LW, 1.5);
    hold on
    semilogy(tt(200:400:end), EigDBO(:,200:400:end),'x','color', RGB2, LW, 1.5);
end

h= zeros(5,1);
h(1) = plot(NaN, NaN, '--k',LW, 1.5);
h(2) = plot(NaN, NaN, ':','color', RGB1, LW, 2.5);
h(3) = plot(NaN, NaN, '-d','color', RGB3, LW, 1.5,'MarkerSize',3.5);
h(4) = plot(NaN, NaN, '-.d','color', RGB5, LW, 1.5,'MarkerSize',3.5);
h(5) = plot(NaN, NaN, '--x','color', RGB2, LW, 1.5);
legend(h,'Analytical','BO','DO', 'PI-DO', 'DBO');
xlim([t0 tf])
xlabel('Time')
xticks([0 pi/2 pi 3*pi/2 2*pi 5*pi/2 3*pi])
xticklabels({'0','\pi/2', '\pi', '3\pi/2', '2\pi', '5\pi/2', '3\pi'})
ylabel('Eigenvalues')
set(gca, 'FontSize', 16, 'Fontname', 'Times New Roman');
ylim([10^-12 10^2])
fname = sprintf('ErrorPlots/%s/Eigenvalues',epname);
saveas(gcf,fname,'epsc')
saveas(gcf, fname,'fig')
figure(5)
if(BO)
    semilogy(tt,Error, ':','color', RGB1, LW, 2.5);
    hold on
end
if(DO)
    semilogy(tt(1:400:end), Errord(1:400:end),'d','color', RGB3,'MarkerSize',3.5);
    hold on
    semilogy(tt, Errord,'-','color', RGB3, LW, 1.5);
end
if(PIDO)
    semilogy(tt(400:400:end), Errordpi(400:400:end),'d','color', RGB5,'MarkerSize',3.5);
    hold on
    semilogy(tt, Errordpi,'-.','color', RGB5, LW, 1.5);
end
if(DBO)
    semilogy(tt, Errorp, '-','color', RGB2, LW, 1.5);
    hold on
    semilogy(tt(100:400:end), Errorp(100:400:end),'x','color', RGB2, LW, 1.5);
end
h= zeros(4,1);
h(1) = plot(NaN, NaN, ':','color', RGB1, LW, 2.5);
h(2) = plot(NaN, NaN, '-d','color', RGB3, LW, 1.5,'MarkerSize',3.5);
h(3) = plot(NaN, NaN, '-.d','color', RGB5, LW, 1.5,'MarkerSize',3.5);
h(4) = plot(NaN, NaN, '-x','color', RGB2, LW, 1.5);
legend(h, 'BO','DO','PI-DO', 'DBO');
title(['Mean error'])
xlim([t0 tf])
ylim([10^-12 10^-4])
xlabel('Time')
xticks([0 pi/2 pi 3*pi/2 2*pi 5*pi/2 3*pi])
xticklabels({'0','\pi/2', '\pi', '3\pi/2', '2\pi', '5\pi/2', '3\pi'})
ylabel('$\mathrm{L_2}$ error','interpreter', 'latex')
set(gca, 'FontSize', 16, 'Fontname', 'Times New Roman');
fname = sprintf('ErrorPlots/%s/MeanError',epname);
saveas(gcf,fname,'epsc')
saveas(gcf,fname,'fig')

figure(6)
if(BO)
    semilogy(tt,Error_var, ':','color', RGB1, LW, 2.5);
    hold on
end
if(DO)
    semilogy(tt(1:400:end), Error_vard(1:400:end),'d','color', RGB3,'MarkerSize',3.5);
    hold on
    semilogy(tt, Error_vard,'-','color', RGB3, LW, 1.5);
end
if(PIDO)
    semilogy(tt(400:400:end), Error_vardpi(400:400:end),'d','color', RGB5,'MarkerSize',3.5);
    hold on
    semilogy(tt, Error_vardpi,'-.','color', RGB5, LW, 1.5);
end
if(DBO)
    semilogy(tt, Error_varp, '-','color', RGB2, LW, 1.5);
    hold on
    semilogy(tt(100:400:end), Error_varp(100:400:end),'x','color', RGB2, LW, 1.5);
end
h= zeros(4,1);
h(1) = plot(NaN, NaN, ':','color', RGB1, LW, 2.5);
h(2) = plot(NaN, NaN, '-d','color', RGB3, LW, 1.5,'MarkerSize',3.5);
h(3) = plot(NaN, NaN, '-.d','color', RGB5, LW, 1.5,'MarkerSize',3.5);
h(4) = plot(NaN, NaN, '-x','color', RGB2, LW, 1.5);
legend(h,'BO', 'DO', 'PI-DO', 'DBO');
title(['Variance error'])
xlim([t0 tf])
ylim([10^-12 10^-4])
xlabel('Time')
xticks([0 pi/2 pi 3*pi/2 2*pi 5*pi/2 3*pi])
xticklabels({'0','\pi/2', '\pi', '3\pi/2', '2\pi', '5\pi/2', '3\pi'})
ylabel('$\mathrm{L_2}$ error','interpreter', 'latex')
set(gca, 'FontSize', 16, 'Fontname', 'Times New Roman');
fname = sprintf('ErrorPlots/%s/VarError',epname);
saveas(gcf,fname,'epsc')
saveas(gcf,fname,'fig')
