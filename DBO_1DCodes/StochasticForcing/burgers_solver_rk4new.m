function soln = burgers_solver_rk4new(Ns, x, u0,tini, tf, dt, nu, xi)

NStep = ceil((tf-tini) / dt);
u = zeros(Ns,1);
cu = u0;
u(:,1) = cu(:,1);
for j=1:NStep
    t  = dt*(j-1)+tini;
    t1 = t+dt/2;
    t2 = t+dt/2;
    t3 = dt*j+tini;
    
    k1 = rhs1(cu, Ns, nu, xi, t);
    k2 = rhs1(cu+dt*k1/2, Ns, nu, xi, t1);
    k3 = rhs1(cu+dt*k2/2, Ns, nu, xi, t2);
    k4 = rhs1(cu+dt*k3, Ns, nu, xi, t3);
    
    u(:,1) = cu + dt*(k1+2*k2+2*k3+k4)/6;
    
    cu(:,1) = u(:,1);
end

soln = u;

end