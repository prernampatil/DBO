function f=Y_exact(t)

global xr
global wr
global epsilon
% ae(1) = sin(t);
% ae(2) = cos(3*t);

ae(1) = 4.5 + sin(t);
ae(2) = 1.5 + cos(3*t);

% ae(1) = 1.5  + sin(t);
% ae(2) = 1.5 + cos(3*t);

% ae(1) = 1-sin(4*t+6*pi/20);
% ae(2) = 1+cos(t);
% ae(1) = sin(t).^2 +1.0;
% ae(2) = 1.0-1.0*tanh(5.0*t)+1.5*10^-5;
scale = 1.0;

Ye = [];

% Ye = [Ye scale*ae(1)*sin(pi*xr(:,1)-t)*sqrt(2)];
% Ye = [Ye scale*ae(2)*cos(pi*xr(:,2)-t)*sqrt(2)];
Ye = [Ye scale*ae(1)*sin(pi*xr(:,1)-t)*sqrt(2)];
Ye = [Ye scale*ae(2)*(epsilon*cos(pi*xr(:,2)-t))*sqrt(2)];

f=Ye;