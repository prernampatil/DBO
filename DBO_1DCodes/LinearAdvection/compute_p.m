function p = compute_p(Lu, Y, w, dubar, du)
% Compute p(k,j) = p_j(x_k) = E[L[u]Y_j](x_k), k=1,...,Ns, j=1,...,N.
global vmean 
global sigma
global xr
	Ns = size(Lu,1);
	N  = size(Y,2);
	p  = zeros(Ns,N);
    
	for i=1:N
		p_test(:,i) = Lu*(Y(:,i).*w);
        
    end

end


