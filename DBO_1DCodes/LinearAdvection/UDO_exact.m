function f=UDO_exact(t)

global xx; 
global vmean; 

Ue = []; 
Ue = [Ue sin(xx - vmean*t*pi)/sqrt(pi)];
Ue = [Ue -cos(xx - vmean*t*pi)/sqrt(pi)];

f=Ue;