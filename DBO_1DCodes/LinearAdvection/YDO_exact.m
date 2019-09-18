function f=YDO_exact(t)

global xr
global sigma;

Ye =[];
Ye = [Ye sqrt(pi)*(cos(sigma*xr*t*pi) - sin(sigma*t*pi)/(pi*sigma*t))];
Ye = [Ye sqrt(pi)*sin(sigma*xr*t*pi)];

f=Ye;