function C = tree2sectChannel(T,th,th0,ifRnd,ifFlip)
% TREE2SECTCHANNEL - generate a (hydrodynamical) channel in a circular
% sector domain based on the topology of a tree T.
% 
% Input:
%   T: represent a tree with a vector of parent pointers. 
%       (e.g. T = FibTree(5) is a Fibonacci tree of height 5)
%   th: angular range of the target circular sector domain. (The sector is
%       bounded by [-th, th].)
%   th0: optional angular shift, so that the angular range of the sector is
%       [th0 - th, th0 + th].
% Output: 
%   C = [x(1) x(2) ... x(n); y(1) y(2) ... y(n)], coordinates of channel
%       boundary, where components of the boundary are separated by NaN's.
%
% Requires: tree2channel.m, rect2sect.m, FibTree.m
% Bowei Wu, 10-15-2018

if nargin == 0, test(); return; end
if nargin < 5, ifFlip = 0; end
if nargin < 4 || isempty(ifRnd), ifRnd = 1; end
if nargin < 3, th0 = 0; end
if nargin < 2, th = pi/3; end

C = tree2channel(T);
if ifFlip, C(1,:) = 1 - C(1,:); C = fliplr(C); end % mirror reflection of the channel

% transform into a right-hand pointing tree bounded by [0 1, -th, th]
x = 1-C(2,:);
y = C(1,:);
x = (x - min(x))/(max(x)-min(x));
y = 2*th*(y - min(y))/(max(y)-min(y)) - th;

[xout,yout] = rect2sect(x,y,th0,ifRnd);
C = [xout;yout];
end

function test()
T = FibTree(5);
C = tree2sectChannel(T);
plot(C(1,:),C(2,:))
axis equal
axis([-3,3,-3,3])

hold on
r = max(C(1,:))*1.02;
t = linspace(0,1,200);
plot(r*exp(2*pi*1i*t),'--k')
hold off
end
