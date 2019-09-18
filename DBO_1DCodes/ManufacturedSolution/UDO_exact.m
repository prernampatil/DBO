function f=UDO_exact(t)

global s1x
global c1x
global s2x
global c2x

Ue = [];
Ue = [Ue [c1x s1x]*[cos(t) ; sin(t)]/sqrt(pi)];
Ue = [Ue [c2x s2x]*[cos(3*t) ; sin(3*t)]/sqrt(pi)];

f=Ue;