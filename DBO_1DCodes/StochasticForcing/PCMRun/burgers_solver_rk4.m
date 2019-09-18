function soln = burgers_solver_rk4(Ns, x, u0, tf, dt, nu, xi)

NStep = ceil(tf / dt);
u = zeros(Ns,NStep+1);         % solution. u(i,j) = u(t_j, x_i)
cu = reshape(u0, Ns,1);
u(:,1) = cu;
for j=1:NStep
    t  = dt*(j-1);
    t1 = t+dt/2;
    t2 = t+dt/2;
    t3 = dt*j;
    
    k1 = rhs1(cu, Ns, nu, xi, t);
    k2 = rhs1(cu+dt*k1/2, Ns, nu, xi, t1);
    k3 = rhs1(cu+dt*k2/2, Ns, nu, xi, t2);
    k4 = rhs1(cu+dt*k3, Ns, nu, xi, t3);
    
    u(:,j+1) = cu + dt*(k1+2*k2+2*k3+k4)/6;
    
    cu = u(:,j+1);
end

soln = u;

end