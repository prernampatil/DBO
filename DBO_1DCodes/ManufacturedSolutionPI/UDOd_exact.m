function f=UDOd_exact(t)

global s1x
global c1x
global s2x
global c2x

Uet =[];
Uet = [Uet [c1x s1x]*[-sin(t) ; cos(t)]/sqrt(pi)];
Uet = [Uet [c2x s2x]*[-3*sin(3*t) ; 3*cos(3*t)]/sqrt(pi)];

f=Uet;