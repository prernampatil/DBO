% Code to solve the Dirichlet boundary conditions
% Linear advection diffusion equation
% Updated to solve using the DO method as well
% Written by: Prerna Patil
% Date: 24th April 2021
% Last updated: 3rd May 2021

clc
clear
clear global
close all

rng('default')
set(0,'defaulttextinterpreter','latex')
LW = 'linewidth';
RGB1 = [0 113 188 ]/norm([0 113 188 ]);
RGB2 = [216 82 24 ]/norm([216 82 24 ]);
RGB3 = [20 82 24  ]/norm([20  82 24 ]);
global NModes
for NModes = 9
    clearvars -except RGB1 RGB2 RGB3 LW NModes
    %% -------------- Global variables --------------
    global LVel D1 D2 nu np wr wp dt sigmat sigmax
    %% -------------- Time Discretization ---------------
    dt        = 0.0005;
    T0        = 0;
    Ts        = dt;            % Final time
    Tf        = 5;
    nTimeStep = ceil((Tf-T0) / dt);
    nu        = 0.05;
    np        = 4;
    nFE       = 101;
    sigmat    = 1;
    sigmax    = 1;
    %% -------------- Space Discretization ---------------
    Xmin               = 0; Xmax = 5;
    fe_mesh            = linspace(Xmin, Xmax, nFE+1);
    [~, x, we, wp, De] = gen_global_coordinate_system(np, fe_mesh);
    we                 = we(1:np+1);
    N                  = length(x);
    
    %% Setting the boundary condition
    bc_type.left  = 'D';
    bc_type.right = 'N';
    bc = [0 0];
    %% -------------- Random Discretization ---------------
    xi = load('ColPnts_d8_333.dat');
    wr = load('ColWgts_d8_333.dat');
    Nr = size(xi,1);
    Nd = size(xi,2);
    LVel   = 1;
    %% ---------------- Differentiation Matrices -----------------
    bc_left = bc_type.left;
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
    
    % Set up the stochastic initial conditions (DBO)
    Y = xi(:,1:Nd);
    U = Ux(:,1:Nd);
    S = diag(Lx(1:Nd));
    
    % Set up the initial conditions for the PCM solver
    Phi = cos(2*pi*x) + sigmax*U*S*Y';
    % Enforce the boundary conditions for PCM
    Phi(1,:) = gkern(xi,Ut,L,1);
    Phi(end,:) = ( -D1(end,end-np:end-1)*Phi(end-np:end-1,:))/D1(end,end);
    
    % Compute the SVD of the solution
    [U,S,Y] = ComputeKL(Phi,wp,wr);
    
    % Obtain the DO solution
    UDO = U;
    YDO = Y*S;
    
    % Set parameters for error computation
    nSkip = 1;
    ErrDBO     = zeros(1,nTimeStep);
    ErrDO      = zeros(1,nTimeStep);
    ErrBndry   = zeros(1,nTimeStep);
    ErrBndryDO = zeros(1,nTimeStep);
    Sigma   = zeros(NModes,nTimeStep);
    SDO     = zeros(NModes,nTimeStep);
    SDNS    = zeros(NModes,nTimeStep/nSkip);
    
    % Save the values for the BC of the modes
    BCModesKL  = zeros(NModes,nTimeStep);
    BCModesDBO = zeros(NModes,nTimeStep);
    BCModesDO  = zeros(NModes,nTimeStep);
    
    for m = 1:nTimeStep-4
        % RK4 time integration
        t = (m-1)*dt;
        T(m) = t+dt;
        gdot     = gkerndot(xi,Ut,L,m);
        k1_Phi   = fPhi(t, Phi, gdot);
        [k1_U, k1_S, k1_Y] = fDBO(t, U, S, Y, gdot);
        [k1_UDO, k1_YDO  ] = fDO(t,UDO, YDO,gdot);
        
        gdot     = gkerndot(xi,Ut,L,m);
        k2_Phi   = fPhi(t+0.5*dt, Phi+0.5*dt*k1_Phi,gdot);
        [k2_U, k2_S, k2_Y] = fDBO(t+0.5*dt, U+0.5*dt*k1_U, S+0.5*dt*k1_S, Y+0.5*dt*k1_Y, gdot);
        [k2_UDO, k2_YDO  ] = fDO(t+0.5*dt ,UDO+0.5*dt*k1_UDO, YDO+0.5*dt*k1_YDO,gdot);
        
        gdot     = gkerndot(xi,Ut,L,m);
        k3_Phi   = fPhi(t+0.5*dt, Phi+0.5*dt*k2_Phi,gdot);
        [k3_U, k3_S, k3_Y] = fDBO(t+0.5*dt, U+0.5*dt*k2_U, S+0.5*dt*k2_S, Y+0.5*dt*k2_Y, gdot);
        [k3_UDO, k3_YDO  ] = fDO(t+0.5*dt ,UDO+0.5*dt*k2_UDO, YDO+0.5*dt*k2_YDO,gdot);
        
        gdot     = gkerndot(xi,Ut,L,m);
        k4_Phi   = fPhi(t+dt,     Phi+dt*k3_Phi,    gdot);
        [k4_U, k4_S, k4_Y] = fDBO(t+dt,     U+dt*k3_U,     S+dt*k3_S,     Y+dt*k3_Y, gdot);
        [k4_UDO, k4_YDO  ] = fDO(t+dt ,UDO+dt*k3_UDO, YDO+dt*k3_YDO,gdot);
        %% ---------%%-----------------
        Phi  = Phi + dt*(k1_Phi + 2.0*k2_Phi + 2.0*k3_Phi + k4_Phi)/6;
        % DBO update
        U    = U   + dt*(k1_U   + 2.0*k2_U   + 2.0*k3_U   + k4_U)/6;
        S    = S   + dt*(k1_S   + 2.0*k2_S   + 2.0*k3_S   + k4_S)/6;
        Y    = Y   + dt*(k1_Y   + 2.0*k2_Y   + 2.0*k3_Y   + k4_Y)/6;
        % DO update
        UDO    = UDO   + dt*(k1_UDO   + 2.0*k2_UDO   + 2.0*k3_UDO   + k4_UDO)/6;
        YDO    = YDO   + dt*(k1_YDO   + 2.0*k2_YDO   + 2.0*k3_YDO   + k4_YDO)/6;
        
        Err = Phi - U*S*Y';
        ErrDBO(m) = sqrt((wp)'*((Err.*Err)*(wr)));
        Err = Phi - UDO*YDO';
        ErrDO(m)  = sqrt((wp)'*(Err.*Err)*(wr));
        
        % Rotation of modes DBO
        [RU,S,RV]  = svd(S);
        Sigma(:,m) = diag(S);
        U = U*RU;
        Y = Y*RV;
        
        % Rotate the DO solution to KL
        C0 = YDO'*diag(wr)*YDO;
        [RUDO, sdo, RVDO ] = svd(C0);
        UDO = UDO*RUDO;
        YDO = YDO*RVDO;
        % Compute error at boundary
        Err = gkern(xi,Ut,L,m) - U(1,:)*S*Y';
        ErrBndry(m) = sqrt(((Err.*Err)*wr));
        Err = gkern(xi,Ut,L,m) - UDO(1,:)*YDO';
        ErrBndryDO(m) = sqrt(((Err.*Err)*wr));
        
        if(mod(m,nSkip)==0)
            [UKL,SKL,YKL]= ComputeKL(Phi,wp,wr);
            SDNS(:,m/nSkip) = diag(SKL);
            SDO(:,m/nSkip)  = sqrt(diag(sdo));
        end
        BCModesKL(:,m)  = UKL(1,:);
        BCModesDBO(:,m) = U(1,:);
        BCModesDO(:,m)  = UDO(1,:);
        
        if(mod(m,200)==0)
            % Plot the evolution of 4 samples
            figure(1)
            plot(x,Phi(:,1:5:20),LW,1.5);
            drawnow
        end
    end
    
    tempstr = sprintf('Dirichlet_d%d_r%d_dt.mat',Nd,NModes);
    save(tempstr,'LW','RGB1','RGB2','RGB3','T','dt','nTimeStep',...
        'ErrDBO','ErrDO','nSkip','SDNS','Sigma','BCModesKL','BCModesDBO','BCModesDO',...
        'ErrBndry','ErrBndryDO','NModes','SDO','Ts');
    
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
SingularValError;

