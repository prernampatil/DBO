function f=ud_exact(t)

global s1x
global c1x
global s2x
global c2x

uet =[s1x c1x]*[-sin(t) ; -cos(t)];

f=uet;