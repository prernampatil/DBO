function [F_U, F_S, F_V] = fDBO(t, U, S, V, gdot)
global wr wp

%% Data driven version 
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
global D1 D2 np LVel nu alpha
F_U = (-LVel*(D1*U*S) + nu*D2*U*S);
% Set the Robin boundary at x=0
F_U(1,:)= (gdot -(D1(1,2:np+1)*F_U(2:np+1,:)*V'))/(D1(1,1)+alpha)*diag(wr)*V;
% Set the Neumann boundary at x=5
F_U(end,:)= (-(D1(end,end-np:end-1)*F_U(end-np:end-1,:)*V'))/D1(end,end)*diag(wr)*V; 

F_V = V*(wp.*F_U)'*U;
F_V = F_V - V* (V'*diag(wr)*F_V);
F_S = (wp.*U)'*F_U;
F_U = (F_U - U* ((wp.*U)'*F_U))/S;

