function gdot = gdot(t,xi)
% 2D random space 
% gdot = 0.01*3*pi*(cos(2*pi*xi(:,1)) + sin(4*pi*xi(:,2)))'*cos(3*pi*t);

% 8D random space 
global nu 
a = 1; 
Nd = size(xi,2);
sigma = 0.01;
gdot = exp(-nu*t)*(-nu*sin(2*pi*a*t) + 2*pi*a*cos(2*pi*a*t));
for i = 0:Nd/2 -1
   gdot = gdot + sigma*2*pi*(2*i+1)*cos(2*pi*(2*i+1)*t)*xi(:,2*i+1)'...
               + sigma*2*pi*(2*i+2)*sin(2*pi*(2*i+2)*t)*xi(:,2*i+2)';
end
