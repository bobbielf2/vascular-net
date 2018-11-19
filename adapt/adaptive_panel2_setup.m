function [p, stdpan] = adaptive_panel2_setup(tol)
% Setup standard panel for use in adaptive split test
% Input:
%  tol - relative desired tolerance
% Outputs:
%  p - num panel nodes: this is not chosen by user, rather, set based on tol.
%  stdpan - struct for standard panel
% Barnett 10/29/18
p = ceil(2 + log10(1/tol));  % guess
[stdpan.nodes, stdpan.weights] = gauss(p);  % col vec nodes on [-1,1]
stdpan.testpts = linspace(-1,1,p)';   % col vec test pts on [-1,1]
stdpan.interpmat = interpmat_1d(stdpan.testpts,stdpan.nodes);
