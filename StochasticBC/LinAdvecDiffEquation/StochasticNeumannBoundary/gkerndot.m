function gdot = gkerndot(xi,Ut,L,m)

global dt sigmat

Utdot = 1/(12*dt)*(-25*Ut(m,:) + 48*Ut(m+1,:) -36*Ut(m+2,:) + 16*Ut(m+3,:)- 3*Ut(m+4,:));
% Utdot = 1/(12*dt)*(-3*Ut(m-1,:) - 10*Ut(m,:) +18*Ut(m+1,:) -6*Ut(m+2,:) + Ut(m+3,:));

t = m*dt;
% Derivative of the Kernel
gdot = (2.0*pi)^2*cos(2*pi*t)+ sigmat*Utdot*diag(L)*xi';
