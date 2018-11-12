function T = FibTree(h)
% FIBTREE - Generate a Fibonacci tree of height h
%
% Input:
%   h: height of the tree
%   T: a vector of parent pointers. The parent of Node k is the Node T(k)
%
% Bowei Wu, 10-15-2018

if nargin == 0, test(); return; end
if h == 0
    T = [0];
elseif h == 1
    T = [0 1];
else
    Th_1 = FibTree(h-1);
    Th_2 = FibTree(h-2);
    T = [0 Th_1+1 2 Th_2(2:end)+length(Th_1)+1];
    % T_h = [0 T_{h-1}+1 2 T_{n-2}(2:end)+c_{n-1}], where c_n = #nodes in T_n
end
end

function test()
T = FibTree(5);
lf = setdiff(1:numel(T),T); % leaf nodes
[x,y] = treelayout(T);
treeplot(T)
hold on
scatter(x(lf),y(lf),'g')
hold off
title('Fibonacci Tree')
end
