function f = fPhi(t, Phi, gdot)
global D1 D2 LVel nu np alpha
f = -LVel*D1*Phi + nu*D2*Phi;
% Set the Robin boundary at x=0
f(1,:)= (gdot - D1(1,2:np+1)*f(2:np+1,:))/(alpha+D1(1,1)); 
% Set the Neumann boundary at x=L
f(end,:) = (-D1(end,end-np:end-1)*f(end-np:end-1,:))/D1(end,end);

