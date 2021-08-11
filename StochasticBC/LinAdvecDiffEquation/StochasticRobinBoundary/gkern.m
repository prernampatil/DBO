function g = gkern(xi,Ut,L,m)
global dt sigmat 
t = m*dt;
g = (-cos(2*pi*t)+(2*pi)*sin(2*pi*t)) + sigmat*Ut(m,:)*diag(L)*xi';