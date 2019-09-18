function f=u_exact(t)

global xx
global vmean 
global sigma 

% Solution mean 
ue = sin(xx - pi*vmean*t)*sin(pi*sigma*t)/(pi*sigma*t);

f=ue;