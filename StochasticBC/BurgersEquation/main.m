% Code to solve the burger's equation
% using DO and DBO
% u_t = -u*u_x + nu*u_xx
% Written by: Prerna Patil
% Last updated: 7th may 2021

clc
clear
clear global
close all

rng('default')
set(0,'defaulttextinterpreter','latex')
LW   = 'linewidth';
RGB1 = [0 113 188 ]/norm([0 113 188 ]);
RGB2 = [216 82 24 ]/norm([216 82 24 ]);
RGB3 = [20 82 24  ]/norm([20  82 24 ]);
global NModes
for NModes = 4:2:8
    clearvars -except NModes RGB1 RGB2 RGB3 LW
    %% -------------- Global variables --------------
    global D1 D2 nu np wr wp dt sigmat sigmax
    %% -------------- Time Discretization ---------------
    dt        = 0.00025;
    T0        = 0;
    Ts        = 0.3;            
    Tf        = 5;          % Final time
    nTimeStep = ceil((Tf-T0) / dt);
    nTswitch  = floor((Ts-T0) / dt);
    nu        = 0.05;
    np        = 4;
    nFE       = 101;
    sigmax    = 0.005;
    sigmat    = 0.01;
    
    %% -------------- Space Discretization ---------------
    Xmin               = 0; Xmax = 1;
    fe_mesh            = linspace(Xmin, Xmax, nFE+1);
    [~, x, we, wp, De] = gen_global_coordinate_system(np, fe_mesh);
    we                 = we(1:np+1);
    N                  = length(x);
    
    %% Setting the boundary condition
    bc_type.left  = 'D';
    bc_type.right = 'N';
    bc = [0 0];
    
    %% -------------- Random Discretization ---------------
    xi = load('ColPnts_d4_256.dat'); % Points and weights loaded from the MePCM code
    wr = load('ColWgts_d4_256.dat');
    Nr = size(xi,1);
    Nd = size(xi,2);
    
    
    %% ---------------- Initial Condition for DBO and DO -----------------
    Ydbo = zeros(Nr,NModes);
    Udbo = zeros(N,NModes);
    Sdbo = zeros(NModes,NModes);
    
    Ydo  = zeros(Nr, NModes);
    Udo  = zeros(N,NModes);
    
    %% ---------------- Differentiation Matrices -----------------
    bc_left  = bc_type.left;
    bc_right = bc_type.right;
    D1 = gen_advection_matrix(np, fe_mesh);
    D1 = diag(wp)\D1;
    D2 = -gen_diff_matrix(np, nFE, we, De, nu, bc_left, bc_right); D2=diag(wp)\D2;
    
    %% Load the kernel file
    % Temporal kernel
    % Run the File: GaussianKernelTime.m to generate the *.m file
    tempstr = sprintf('GaussKern_d%d_Tf%d.mat',Nd,Tf);
    load(tempstr);
    L = sqrt(L);
    %% ---------------- Initial Condition -----------------
    % Import the Spatial Gaussian kernel
    % Run the File: GaussianKernelSpatial.m to generate the *.m file
    tempstr = sprintf('SpatialGaussKern_d%d.mat',Nd);
    load(tempstr);
    Lx = sqrt(Lx);

    % Set up the stochastic initial conditions
    Y = xi(:,1:Nd);
    U = Ux(:,1:Nd);
    S = diag(Lx(1:Nd));
    
    % Set the initial conditions for the PCM solver
    u        = sin(2*pi*x) + sigmax*U*S*Y';
    
    % Enforce the boundary conditions for the PCM solver
    u(1,:)   = (gkern(xi,Ut,L,1));
    u(end,:) = (-D1(end,end-np:end-1)*u(end-np:end-1,:))/D1(end,end);
    
    %% ----------------- Run the code till stochasticity develops -----------------
    for m = 1:nTswitch
        % RK4 time integration
        t        = (m-1)*dt;
        T(m)     = t+dt;
        gdot     = gkerndot(xi,Ut,L,m);
        k1_Phi   = fPhi(t, u, gdot);
        
        gdot     = gkerndot(xi,Ut,L,m);
        k2_Phi   = fPhi(t+0.5*dt, u+0.5*dt*k1_Phi, gdot);
        
        gdot     = gkerndot(xi,Ut,L,m);
        k3_Phi   = fPhi(t+0.5*dt, u+0.5*dt*k2_Phi, gdot);
        
        gdot     = gkerndot(xi,Ut,L,m+1);
        k4_Phi   = fPhi(t+dt,     u+dt*k3_Phi, gdot);
        
        u      = u + dt*(k1_Phi + 2.0*k2_Phi + 2.0*k3_Phi + k4_Phi)/6;
        
    end
    
    % Compute the SVD of the solution
    [Udbo,Sdbo,Ydbo] = ComputeKL(u,wp,wr);
    % Get IC for DO
    UDO = Udbo;
    YDO = Ydbo*Sdbo';
    % Set parameters for error computation
    nSkip = 1;
    ErrDBO     = zeros(1,nTimeStep-nTswitch);
    ErrDO      = zeros(1,nTimeStep-nTswitch);
    ErrBndry   = zeros(1,nTimeStep-nTswitch);
    ErrBndryDO = zeros(1,nTimeStep-nTswitch);
    Sigma      = zeros(NModes,nTimeStep);
    SDNS       = zeros(NModes,nTimeStep/nSkip);
    SDO        = zeros(NModes,nTimeStep/nSkip);
    % Save the values for the BC of the modes
    BCModesKL  = zeros(NModes,nTimeStep-2);
    BCModesDBO = zeros(NModes,nTimeStep-2);
    BCModesDO  = zeros(NModes,nTimeStep-2);
    
    for m = nTswitch+1:nTimeStep-3
        % RK4 time integration
        t = (m-1)*dt;
        T(m) = t+dt;
        gdot     = gkerndot(xi,Ut,L,m);
        k1_Phi   = fPhi(t, u, gdot);
        [k1_Udbo, k1_Sdbo, k1_Ydbo] = fDBO(t, Udbo, Sdbo, Ydbo, gdot);
        [k1_UDO, k1_YDO] = fDO(t, UDO, YDO, gdot);
        
        gdot     = gkerndot(xi,Ut,L,m);
        k2_Phi   = fPhi(t+0.5*dt, u+0.5*dt*k1_Phi,gdot);
        [k2_Udbo, k2_Sdbo, k2_Ydbo] = fDBO(t+0.5*dt, Udbo+0.5*dt*k1_Udbo, Sdbo+0.5*dt*k1_Sdbo, Ydbo+0.5*dt*k1_Ydbo, gdot);
        [k2_UDO, k2_YDO] = fDO(t+0.5*dt, UDO+0.5*dt*k1_UDO, YDO+0.5*dt*k1_YDO, gdot);
        
        gdot     = gkerndot(xi,Ut,L,m);
        k3_Phi   = fPhi(t+0.5*dt, u+0.5*dt*k2_Phi,gdot);
        [k3_Udbo, k3_Sdbo, k3_Ydbo] = fDBO(t+0.5*dt, Udbo+0.5*dt*k2_Udbo, Sdbo+0.5*dt*k2_Sdbo, Ydbo+0.5*dt*k2_Ydbo, gdot);
        [k3_UDO, k3_YDO] = fDO(t+0.5*dt, UDO+0.5*dt*k2_UDO, YDO+0.5*dt*k2_YDO, gdot);
        
        gdot     = gkerndot(xi,Ut,L,m);
        k4_Phi   = fPhi(t+dt,     u+dt*k3_Phi,    gdot);
        [k4_Udbo, k4_Sdbo, k4_Ydbo] = fDBO(t+dt,     Udbo+dt*k3_Udbo,     Sdbo+dt*k3_Sdbo,     Ydbo+dt*k3_Ydbo, gdot);
        [k4_UDO, k4_YDO] = fDO(t+dt, UDO+dt*k3_UDO, YDO+dt*k3_YDO, gdot);
        
        %% ---------%%-----------------
        u  = u + dt*(k1_Phi + 2.0*k2_Phi + 2.0*k3_Phi + k4_Phi)/6;
        % DBO update
        Udbo    = Udbo   + dt*(k1_Udbo   + 2.0*k2_Udbo   + 2.0*k3_Udbo   + k4_Udbo)/6;
        Sdbo    = Sdbo   + dt*(k1_Sdbo   + 2.0*k2_Sdbo   + 2.0*k3_Sdbo   + k4_Sdbo)/6;
        Ydbo    = Ydbo   + dt*(k1_Ydbo   + 2.0*k2_Ydbo   + 2.0*k3_Ydbo   + k4_Ydbo)/6;
        % DO update
        UDO    = UDO   + dt*(k1_UDO   + 2.0*k2_UDO   + 2.0*k3_UDO   + k4_UDO)/6;
        YDO    = YDO   + dt*(k1_YDO   + 2.0*k2_YDO   + 2.0*k3_YDO   + k4_YDO)/6;
        
        Err = u - Udbo*Sdbo*Ydbo';
        ErrDBO(m) = abs((wp)'*(Err.*Err)*(wr));
        Err = u - UDO*YDO';
        ErrDO(m)  = abs((wp)'*(Err.*Err)*(wr));
        
        % Rotation of modes
        [RU,Sdbo,RV]  = svd(Sdbo);
        Sigma(:,m) = diag(Sdbo);
        Udbo = Udbo*RU;
        Ydbo = Ydbo*RV;
        
        % Rotate the DO solution to KL
        C0 = YDO'*diag(wr)*YDO;
        [RUDO, sdo, RVDO ] = svd(C0);
        UDO = UDO*RUDO;
        YDO = YDO*RVDO;
        % Compute error at boundary
        Err = gkern(xi,Ut,L,m) - Udbo(1,:)*Sdbo*Ydbo';
        ErrBndry(m) = sqrt(abs((Err.*Err)*wr));
        Err = gkern(xi,Ut,L,m) - UDO(1,:)*YDO';
        ErrBndryDO(m) = sqrt(abs((Err.*Err)*wr));
        
        if(mod(m,nSkip)==0)
            [UKL,SKL,YKL]= ComputeKL(u,wp,wr);
            SDNS(:,m/nSkip) = diag(SKL);
            SDO(:,m/nSkip)  = sqrt(diag(sdo));
        end
        BCModesKL(:,m)  = UKL(1,:);
        BCModesDBO(:,m) = Udbo(1,:);
        BCModesDO(:,m)  = UDO(1,:);
        
        if(mod(m,200)==0)
            % Plot the evolution
            figure(1)
            plot(x,u,LW,1.5);
            drawnow
        end
    end
    
    tempstr = sprintf('BurgersDBC_d%d_r%d.mat',Nd,NModes);
    save(tempstr,'LW','RGB1','RGB2','RGB3','T','dt','nTswitch','nTimeStep',...
        'ErrDBO','ErrDO','nSkip','SDNS','Sigma','BCModesKL','BCModesDBO','BCModesDO',...
        'ErrBndry','ErrBndryDO','NModes','SDO','Ts','Tf');
end

% Post processing of the results
% Plot the in the singular values
plotSingularValues;

% Plot the Global Error
plotGlobalErrors;

% Plot the boundary Error
plotBoundaryErrors;

% Plot evolution of the Boundary modes 
plotBCModes;

% Plot the error in the singular values 
% SingularValError;