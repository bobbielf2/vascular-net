function z = ppfun(pp,t)
% evaluate piecewise polynomial, reparameterized by t in [0,2pi]
%
% Input:
%   pp      piecewise polynomial, see also: mkpp, ppval
%   t       query points, assume t in [0,2pi]

if nargin == 0, test_ppfun; return; end

a = pp.breaks(1);
b = pp.breaks(end);
% t = linspace(0,2*pi,1000);
x = (b-a)/(2*pi)*t+a;
z = ppval(pp,x);


function test_ppfun

breaks=[0,1,2,3];
coefs = [1,-2,2
        1,-1,1
        -1,0.2,1];
pp=mkpp(breaks,coefs);
t = linspace(0,2*pi,1000);
z = ppfun(pp,t);
plot(t,z)