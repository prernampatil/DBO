function [F_U, F_S, F_V] = fDBO(t, U, S, V, gdot)
% global wr wp 

% %% Data driven version 
% F = fPhi(t,U*S*V',xi);
% 
% F_U = F*diag(wr)*V/S;
% F_U = F_U - U*((wp.*U)'*F_U);
% 
% F_V = (wp.*F)'*U/S';
% F_V = F_V - V* (V'*diag(wr)*F_V);
% 
% F_S = (wp.*U)'*F*diag(wr)*V;


%% Model driven version (compressed form) %%
% global D1 D2 np LVel nu
% F_U = (-LVel*(D1*U*S) + nu*D2*U*S);
% % Set the Neumann boundary at x=0
% F_U(1,:)= gdot*diag(wr)*V; 
% % Set the Neumann boundary at x=0
% F_U(end,:)= (-(D1(end,end-np:end-1)*F_U(end-np:end-1,:)*V'))/D1(end,end)*diag(wr)*V; 
% 
% 
% F_V = V*(wp.*F_U)'*U;
% F_V = F_V - V* (V'*diag(wr)*F_V);
% F_S = (wp.*U)'*F_U;
% F_U = (F_U - U* ((wp.*U)'*F_U))/S;



global  wp wr D1 D2 nu np

%% --------- Brute-Force Approach (akin to Data-Driven) -----------
% F = f(t,U*S*V');
%% --------- Model-Driven Approach -----------
% F = -(U*S*V').*(D1*U*S*V') + D2*U*S*V';
% With the mean separated part 
F =   -(U*S*V').*(D1*U*S*V') + nu*D2*U*S*V';
F(1,:) = gdot; %(Enforcing the boundary condition)
% Set the Neumann boundary at x=L
F(end,:)= (-(D1(end,end-np:end-1)*F(end-np:end-1,:)))/D1(end,end);

F_S = U'*diag(wp)*F*diag(wr)*V;
F_U = (F*diag(wr)*V -U*F_S)/S;

F_V = S\(U'*diag(wp)*F-F_S*V');
F_V = F_V';
