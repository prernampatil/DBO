function f=u_exact(t)

global s1x
global c1x
global s2x
global c2x

ue = [s1x c1x]*[ cos(t) ; -sin(t)];

f=ue;