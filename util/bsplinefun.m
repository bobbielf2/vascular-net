function z = bsplinefun(n,t,P,w,u)
% evaluate weighted/unweighted B-splines function, reparameterized by u in [0,2pi]
%
% Input:
%   n       order of B-splines
%   t       knot vector
%   P       control points
%   w       weights associated to P
%   u       query points, assume u in [0,2pi]

if 1
    z = bspline_wdeboor(n,t,P,w,u/(2*pi)); % weighted (NURBS)
else
    z = bspline_deboor(n,t,P,u/(2*pi)); % unweighted (B-splines)
end
z = z(1,:) + 1i*z(2,:);