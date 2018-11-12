function f = randfun(x,m)
% random smooth function
% m = #Fourior modes
% Ref: Trefethen, Exploring ODEs, SIAM 2018

if nargin == 0, test(); return; end
if nargin < 2, m = 5; end
rng(1);
a = randn(m,1)/sqrt(2*m+1);
b = randn(m,1)/sqrt(2*m+1);

f = zeros(size(x));
for k = 1:m
    f = f + a(k)*cos(2*pi*k*x) + b(k)*sin(2*pi*k*x);
end

f = f/sqrt(2);
end

function test()
x = linspace(0,1,1000);
plot(x,randfun(x,20))
end