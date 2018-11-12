function Z = mobius(Z,a,R)
% Mobius transformation mapping the disk (radius=R) to itself, z=0 is 
% mapped to z=a*R.
%
% Input:
%   C = [x(1),x(2),... ; y(1),y(2),...] coordinates
%   a = any complex number, |a|<1
%   R = radius of disk

if nargin ==0, test(); return; end

Z = (Z + a*R)./(1 + a'/R*Z);
end

function test()
R = 3;
r = linspace(0.3, R,30);
th = linspace(0,2*pi,60);
[r, th] = meshgrid(r,th);
Z = r.*exp(1i*th);
mesh(real(Z) - R - 1,imag(Z),0*Z)

a = 0.2 - 0.5i;
Z = mobius(Z,a,R);
hold on
mesh(real(Z) + R + 1,imag(Z),0*Z)
hold off
axis([-2.5,2.5,-1,1]*R)
axis equal off
end