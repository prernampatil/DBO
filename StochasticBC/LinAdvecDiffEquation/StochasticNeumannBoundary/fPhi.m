function f = fPhi(t, Phi, gkerndot)
global D1 D2 LVel nu np 
f = -LVel*D1*Phi + nu*D2*Phi;
% Set the Neumann boundary at x=0
f(1,:)= (gkerndot -D1(1,2:np+1)*f(2:np+1,:))/D1(1,1); 
% Set the Neumann boundary at x=L (no stochasticity)
f(end,:) = ( -D1(end,end-np:end-1)*f(end-np:end-1,:))/D1(end,end);
