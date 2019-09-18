%
% du(i, j) = du_j/dx(t,x_i)
%

function ret = Diff_rk4(u, order)
global vmean 
global x 
    [m n] = size(u);
    ret = zeros(size(u));
    for i=1:n
        ret(:,i) = fourdifft(  u(:,i), order );
    end

end

