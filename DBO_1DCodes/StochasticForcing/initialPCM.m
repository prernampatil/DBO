function [u_pcm, mean, var]= initialPCM(Ns, Nr, xr, x, u0,tini, dt, mu, wr,ts)
mean = zeros(Ns, 1);
var  = zeros(Ns,1);
for i=1:Nr
    
    %  xi = (xr(i)+1)/2;
    u_pcm(:,i) = burgers_solver_rk4new(Ns, x, u0, tini,ts, 0.1*dt, mu, xr(i));
    mean = mean + u_pcm(:,i) * wr(i);
    var = var + u_pcm(:,i).^2 * wr(i);
    
end
