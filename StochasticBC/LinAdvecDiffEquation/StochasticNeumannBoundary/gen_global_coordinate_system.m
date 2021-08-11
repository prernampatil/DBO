function [x, xg, w, wg, D] = gen_global_coordinate_system(np, fe_mesh)
% This function returs a vector including the coordinates of all the collocation 
% point in the elements defined by the vector fe_mesh 
% 
%   Input       P     -> spectral element order
%             fe_mesh -> vector containing the coordinates of 
%                        the element boundaries 
%
%   Output      x     -> coordinates of all nodes points
%                        (with repeated boundaries)
%               xg    -> global coordinates
%               w     -> quadrature weights
%               wg    -> quadrature weights in global coordinates
%               D     -> element differential matrix




K = length(fe_mesh) - 1; % Number of finite elements
P1 = np + 1; % Number of quadrature points in each element
x = zeros(K*P1, 1); w = x; xg = zeros(K*np+1, 1); wg = xg;

for i = 1:K
    [x1, w1, D] = get_GLL_points_and_D_matrix(np, fe_mesh(i), fe_mesh(i+1));
    x((i-1)*P1+1:i*P1) = x1;      w((i-1)*P1+1:i*P1) = w1;
    xg((i-1)*np+1:i*np+1) = x1;
    wg((i-1)*np+1:i*np+1) = wg((i-1)*np+1:i*np+1) + w1;
end



end