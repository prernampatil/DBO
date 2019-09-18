function f=ud_exact(t,xi)

global vmean;
global sigma; 
global x;

ue = sin(x -(vmean+sigma*xi)*t*pi);
f=ue;