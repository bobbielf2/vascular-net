function P = traverse(T,k)
% TRAVERSE - generate a path visiting all the nodes of a subtree (with Node
% k being its root) of a given tree T.
%
% Input:
%   T: represent a tree with a vector of parent pointers. 
%       (e.g. T = FibTree(5) is a Fibonacci tree of height 5)
%   k: a node of the tree T
% Output: 
%   P = [v(1) v(2) ... v(n)]. Ordered list of nodes of T that visits all
%   the nodes of the subtree with Node k being its root. Whenever a leaf of
%   the tree is visited, a NaN is inserted as a marker.
%
% Requires: FibTree.m
% Bowei Wu, 10-15-2018

if nargin == 0, test(); return; end
P = [];
children = find(T == k);
if numel(children) == 0
    P = [k NaN k]; % if Node k is a leaf, use a NaN to separate the components
    return
end
for c = children
    P = [P k traverse(T,c)];
end
P = [P k];
end

function test()
T = FibTree(2);
treeplot(T)
P = traverse(T,1);
title(num2str(P))
end
