function g = g(t,xi)
% 2D random space 
% g = 0.01*(cos(2*pi*xi(:,1)) + sin(4*pi*xi(:,2)))'*sin(3*pi*t);

% 8D random space 
global nu 
a = 1; 
Nd = size(xi,2);
sigma = 0.01;
g = exp(-nu*t)*sin(2*pi*a*t);
for i = 0:Nd/2 -1
   g = g + sigma*sin(2*pi*(2*i+1)*t)*xi(:,2*i+1)'...
         - sigma*cos(2*pi*(2*i+2)*t)*xi(:,2*i+2)';
end
