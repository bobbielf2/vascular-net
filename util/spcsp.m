function pp = spcsp(x,Y,p,varargin)
% SMOOTHING CUBIC SPLINE WITH PERIODIC CONDITIONS
%
% SYNTAX:
%   pp = spcsp(x,Y,p);
%   pp = spcsp(x,Y,p,w);
%
% spcsp(x,Y,p,w) returns the ppform of a cubic smoothing spline for the
% given data x,Y, with periodic conditions. The smoothing spline 
% approximates, at the data site x(j), the given data value Y(:,j), where 
% j=1:length(x). The data values may be scalars, vectors, matrices, 
% or even ND-arrays. If weigth vector w is not supplied, default weigths to
% 1 are used.
%
% INPUTS:
% x = row vector of breaks (parametric variable of the d-dimensional curve)
% Y = d-by-n matrix of the values to approximate
%       d: dimension of the space
%       n: number of points
% p = smoothing parameter
% w = row vector of weights
%
% The smoothing spline f minimizes
%
%    p * sum_j w(j)^(-2) |Y(:,j) - f(x(j))|^2  +  (1-p) * integral |D^2 f|^2
% 
% where the sum is over j=1:length(X).
%
% OUTPUTS:
% pp = ppform of the approximating spline.
%
%
%Copyrights: Massimo Zanetti 2016.
% https://www.mathworks.com/matlabcentral/fileexchange/59463-smoothing-cubic-splines-with-periodic-conditions
% ################# ALGORITHM #############################################

if nargin == 0, spcsp_test; return; end

x=x(:)';
n = numel(x); % number of breaks
[dim,n2] = size(Y); % space dimension and number of points
Y = Y(:,1:end-1);
%check inputs
if ~isempty(varargin)
    w=(varargin{1}(:))';
else
    w=ones(size(x));
end
if n~=n2
    error('size of x and Y do not match');
end
if n~=numel(w)
    error('numel(x)~=numel(w), check dimensions');
end
% h array of sub-intervals between breaks
h = x(2:end)-x(1:end-1); % n-1 array
% S matrix
aux = [ [h(1:end-1), 0] ; 2*(h+[h(end), h(1:end-1)]) ; [0, h(1:end-1)] ]';
S = spdiags( aux , [-1 0 1] , n-1 , n-1 );
S(1,n-1) = h(end);
S(n-1,1) = h(end);
% V matrix
aux = [ [1./h(1:end-1), 0] ; -1./h-1./[h(end), h(1:end-1)] ; [0, 1./h(1:end-1)] ]';
V = spdiags( aux , [-1 0 1] , n-1 , n-1 );
V(1,n-1) = 1/h(end);
V(n-1,1) = 1/h(end);
% generate weight matrix from w, matrix W:=W^-2
W = spdiags( w(1:end-1)'.^(-2) , 0 , n-1 , n-1 );
% constructing the quadratic form: F(a) = 1/2 a'Ua - v'a
AUX = S\V; %AUX = S^-1 V
U = 2*p*W + 12*(1-p)*V'*AUX;
v = 2*p*W*Y';
h = h';
% minimizing and computing coefficients
a = U\v;
c = 3*AUX*a;
d = ([c(2:end,:);c(1,:)]-c)./(3*h);
b = ([a(2:end,:);a(1,:)]-a)./h - c.*h - d.*(h.^2);
coefs = cat(3,d',c',b',a');
% build spline in pp form: breaks=x
pp = mkpp(x,coefs,dim);

function spcsp_test
clc; clear; close all;
%Run this script to see smoothing cubic splines with periodic conditions
%used to approximate closed curves in 2D and 3D space.
%2D curve approximation: smoothing a noisy ellipse
a = 5; b = 2;
t = 0:.2:2*pi;
%noisy ellipse points
Y = [a*cos(t); b*sin(t)] + rand(2,length(t))-.5;
x = 1:numel(t);
xx = linspace(x(1),x(end),10*length(x)+1);
for p=0.1  
    % 
    pp = spcsp(x,Y,p);
    %     
    val = ppval(pp,xx);
    figure('Color','w');
    plot(Y(1,:),Y(2,:),'r--',val(1,:),val(2,:),'b-','LineWidth',5);
    title('spcsp in 2D');
    set(gca,'FontSize',30);
    grid on;
    legend('Noisy curve','Spline approx');
    xlabel('x');
    ylabel('y');
end
% 3D curve approximation: smoothing warped ellipse 
a = 5; b = 2;
t = 0:.2:2*pi;
z = [1:22,10:-1:1];
%noisy ellipse points
Y = [a*cos(t); b*sin(t); z] + rand(3,length(t))-.5;
x = 1:numel(t);
xx = linspace(x(1),x(end),10*length(x)+1);
for p=0.1  
    % 
    pp = spcsp(x,Y,p);
    %     
    val = ppval(pp,xx);
    figure('Color','w');
    plot3(Y(1,:),Y(2,:),Y(3,:),'r--',val(1,:),val(2,:),val(3,:),'b-','LineWidth',5);
    title('spcsp in 3D');
    set(gca,'FontSize',30);
    grid on;
    legend('Noisy curve','Spline approx');
    xlabel('x');
    ylabel('y');
    zlabel('z');
end