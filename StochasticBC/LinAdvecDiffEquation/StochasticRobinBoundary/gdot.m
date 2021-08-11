function gdot = gdot(t,xi)
gdot = 0.01*2*pi*cos(2*pi*t)*(cos(2*pi*xi(:,1)) + sin(4*pi*xi(:,2)))'*cos(2*pi*t);