function [s, W] = make_quadr_net(tol,smoothing)
% MAKE_QUADR_NET: Generate quadratures for the components of the vascular net.
% Input:
%   tol         tolerance for adaptive panel refinement
%   smoothing   0=no smoothing (polygons), 1=smoothing via cubic spline approx
% Output:
%   s   s{i} - quadrature rule for the i-th component of the network
%   W   W{1} - quadr for inner circular wall. 
%       W{2} - quadr for outer circular wall
%       W{i}.r - radius of the circle
%
% B Wu 11-18-18


v = 1; % toggle for plots
if nargin < 2, smoothing = 0; end
if nargin < 1 || isempty(tol), tol = 1e-12; end
[~,pp,ppf] = example_circular_net(smoothing,v);
addpath(genpath('adapt'))

if v, t = linspace(0,2*pi,1000); figure; end

m = numel(ppf);
s = cell(m-2,1);
W = cell(2,1);
for ip = 1:m
    f = ppf{ip};
    if smoothing <= 1 && ip < m-1
        breaks = pp{ip}.breaks;
        breaks = (breaks - breaks(1))/(breaks(end)-breaks(1))*2*pi;
    end
    if smoothing == 0 && ip < numel(ppf)-1
        corner.tc = breaks;
        corner.lam = 3*ones(size(corner.tc));
    else
        corner = [];
    end
    
    [p, stdpan] = adaptive_panel2_setup(tol);
    paramdomain = [0,2*pi];
    tpan = adaptive_panel2(f,p,paramdomain,tol,corner,stdpan);
    if ip < m-1
        s{ip}.tpan = tpan;
        s{ip}.Z = f;
        s{ip} = quadr(s{ip}, [], 'p');
    else
        ipw = ip + 2 - m;
        W{ipw}.tpan = tpan;
        W{ipw}.Z = f;
        W{ipw} = quadr(W{ipw}, [], 'p');
    end
    if v
        plot(f(t))
        hold on
        plot(f(tpan),'.b','markersize',5)
        axis equal
        axis([-1,1,-1,1]*3)
        drawnow
    end
end
W{1}.r = mean(abs(W{1}.x - mean(W{1}.x)));
W{2}.r = mean(abs(W{2}.x - mean(W{2}.x)));
