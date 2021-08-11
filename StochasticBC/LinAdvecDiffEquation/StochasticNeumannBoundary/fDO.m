function [F_U, F_V] = fDO(t, U, V, gdot)
global wr wp 

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
global D1 D2 np LVel nu
C = V'*diag(wr)*V;
F_U = (-LVel*(D1*U*C) + nu*D2*U*C);
% Set the Neumann boundary at x=0
F_U(1,:)= (gdot -(D1(1,2:np+1)*(F_U(2:np+1,:)/C)*V'))/D1(1,1)*diag(wr)*V; 
% Set the Neumann boundary at x=L
F_U(end,:)= (-(D1(end,end-np:end-1)*F_U(end-np:end-1,:)))/D1(end,end); 

F_V = V*((wp.*(-LVel*(D1*U) + nu*D2*U))'*U);
% F_V = F'*diag(wp)*U;
% F_V = (wp.*F_U)'*U;
% F_U = F*diag(wr)*V;
F_U = (F_U - U* ((wp.*U)'*F_U))/C;
