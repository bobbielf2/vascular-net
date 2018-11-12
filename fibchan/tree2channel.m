function C = tree2channel(T)
% TREE2CHANNEL - generate a (hydrodynamical) channel base on the topology
% of a tree (graph).
% 
% Input:
%   T: represent a tree with a vector of parent pointers. 
%       (e.g. T = FibTree(5) is a Fibonacci tree of height 5)
% Output: 
%   C = [x(1) x(2) ... x(n); y(1) y(2) ... y(n)], coordinates of channel
%       boundary, where components of the boundary are separated by NaN's.
%
% Requires: traverse.m, FibTree.m
% Bowei Wu, 10-15-2018

if nargin == 0, test(); return; end
P = traverse(T,1); % path visiting all nodes of T
[x,y] = treelayout(T); % coordinates of the nodes
% TEST
if 1
    % tree width exponentially decaying with tree height
    w = (y-min(y))/(max(y)-min(y)); % normalize height
    L = 3.5;
    w = exp(-abs(w-1).^0.818*L)*0.25; % Adjust width of channel to make circular channel look nice
    d = (y(1)-y(2))/4; % Adjust width of bifurcation to make circular channel look nice
else
    % tree width algebraically decaying with tree height
    w = (y-min(y))/(max(y)-min(y)); % normalize height
    L = 0.8;
    w = (L*w+1-L).^2*0.12; % Adjust width of channel to make circular channel look nice
    d = (y(1)-y(2))/7; % Adjust width of bifurcation to make circular channel look nice
end
ls = [x - w; y]; % left-shift of coords
rs = [x + w; y]; % right-shift of coords
ds = [x; y - d]; % down-shift of coords
C = [P; P];

% preprocessing
for k = 1:numel(T)
    ind = find(P == k); % ind = the index where Node k appears in P
    counts = numel(ind);
    switch counts
        case 2
            C(1:2,ind(1)) = ls(1:2,k);
            C(1:2,ind(2)) = rs(1:2,k);
        case 3
            C(1:2,ind(1)) = ls(1:2,k);
            C(1:2,ind(2)) = ds(1:2,k);
            C(1:2,ind(3)) = rs(1:2,k);
    end
end
end

function test()
T = FibTree(5);
C = tree2channel(T);
treeplot(T)
hold on
plot(C(1,:),C(2,:))
hold off
title('Fibonacci Channel')
end
