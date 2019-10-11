function f=Yd_exact(t)

global xr
global epsilon

ad_e(1) = cos(t);
ad_e(2) = -3*sin(3*t);

ae(1) = 4.5 + sin(t);
ae(2) = 1.5 + cos(3*t);
scale = 1.0;


Yed = [];

Yed = [Yed scale*ad_e(1)*sin(pi*xr(:,1)-t)*sqrt(2) - scale*ae(1)*cos(pi*xr(:,1)-t)*sqrt(2)];
Yed = [Yed scale*ad_e(2)*epsilon*cos(pi*xr(:,2)-t)*sqrt(2) + scale*ae(2)*epsilon*sin(pi*xr(:,2)-t)*sqrt(2)];

f=Yed;