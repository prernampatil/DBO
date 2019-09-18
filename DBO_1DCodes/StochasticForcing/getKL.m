% MS(:, j) = u(:,t; \xi_j) from PCM solution
% mean = \bar{u}(:,t)
%
% Construct KL decomposition from the above solution and return corresponding u
% and Y for DO decomposition.
%
function [u0, Y0] = getKL(MS, mean, N, wr, wp)

[Ns, M] = size(MS);
C0 = zeros(Ns,Ns);

for i=1:M
    MS(:,i) = MS(:,i) - mean;
end

for i=1:Ns
    for j=1:Ns
        C0(i,j) = sum(MS(i,:) .* MS(j,:) .* wr') * wp(j);
    end
end

[V0,~] = eig(C0);

for m=1:N
    u0(:,m) = V0(:,end-m+1)/sqrt(sum(V0(:,end-m+1).^2.*wp));
    
    for n=1:M
        Y0(n,m) = sum(MS(:,n) .* u0(:,m) .* wp);
    end
end

end
