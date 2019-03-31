function [C, pp, ppf] = example_circular_net(smoothing,v)
%% generate a example (hydrodynamical) channel in a circular domain based on
% the topology of a "Fibonacci tree" (a fractal tree).
% Input:
%   smoothing   toggle: 
%                   0 = no smoothing, 
%                   1 = smoothed periodic cubic spline 
%                   2 = smoothed via partition of unity (test)
%   v       toggle for plots
% Output:
%   C       vertices of polygons, C(1,:) and C(2,:) are x- and y-coordinates,
%           different polygons are separated by NaN's
%   pp      piecewise polynomials approximating the polygons C, 
%           pp{k} represents the k-th boundary component.
%           See also: mkpp, ppval
%   ppf     ppf{k} is the inline funtion evaluating pp{k} on the
%           reparameterized domain [0, 2pi]
%
% Requires: FibTree.m, tree2sectChannel.m, fillout.m, (spcsp.m if smoothing)
% Bowei Wu, 10-15-2018, add piecewise polynomial output 11-2

if nargin < 2, v = 0; end % toggle for plots

setup_vasnet();

if nargin == 0, smoothing = 1; end
H = [5,6,7,6,5]; % num of layers for each main branch
% H = [3,3,3];
T = FibTree(H(end));
m = numel(H);
th0 = pi/m;
C = tree2sectChannel(T,th0-0.01); %plot(C(1,:),C(2,:),'-',LW,1.2,MS,2); hold on
for i = 1:m-1
    T = FibTree(H(i));
    ifFlip = mod(i,1); % add variety by mirror reflection of the tree (test)
    C0 = tree2sectChannel(T,th0-0.01,th0*2*i,[],ifFlip); 
    C = [C C0];
end

% post-processing: connecting the components
ind = find(isnan(C(1,:)),1);
C = [C(:,ind:end),C(:,1:ind)];
ind = find(isnan(C(1,:)));
C1 = [];
for i = 1:numel(ind)-1
    C1 = [C1, C(:,ind(i):ind(i+1)-1),C(:,ind(i)+1)];
end
C = [C1,[nan;nan]];

% Shift the center via Mobius transformation
a = (0.14 - 0.1i)*0.999; R = exp(1);
Z = C(1,:) + 1i*C(2,:);
Z = mobius(Z,a,R);
C = [real(Z);imag(Z)];

% plot the channel
if v
    figure
    c = [1,1,1]/2; % color
    t = linspace(0,2*pi,1000); 
end 
ind = find(isnan(C(1,:)));
pp = cell(numel(ind)-1,1);
ppf = cell(numel(ind)-1,1);
for i = 1:numel(ind)-1
    if smoothing == 1 % smoothed cubic splines
        z = C(:,ind(i)+1:ind(i+1)-1); z = z(1,:)+1i*z(2,:);
        z = interp1(z, 1:.5:length(z)); % include midpoints of polygon edges, make interpolant closer to polygon.
        n = length(z);
        pp{i} = spcsp(1:n,z,0.90);
        ppf{i} = @(t) ppfun(pp{i},t);
    elseif smoothing == 2 % partition of unity smoothing
        z = C(:,ind(i)+1:ind(i+1)-2); z = z(1,:)+1i*z(2,:);
        ppf{i} = @(t) smooth_poly(z,t);
    elseif smoothing == 3 % NURBS smoothing
        z = C(:,ind(i)+1:ind(i+1)-2);
        bn = 6;  % order of B-spline (polynomial degree = bn - 1)
        [bt, bP, bw] = polygon_ctrlPts(z, bn);
        ppf{i} = @(t) bsplinefun(bn,bt,bP,bw,t);
    else
        z = C(:,ind(i)+1:ind(i+1)-1); z = z(1,:)+1i*z(2,:);
        n = length(z)-1;
        breaks = 1:n+1;
        coefs = zeros(n,2);
        for k = 1:n
            coefs(k,1) = z(k+1)-z(k);
            coefs(k,2) = z(k);
        end
        pp{i} = mkpp(breaks,coefs);
        ppf{i} = @(t) ppfun(pp{i},t);
    end
    if v
        zz = ppf{i}(t); x = real(zz); y = imag(zz); fill(x,y,c); hold on
    end
end

% internal (source) boundary
r = 0.12;
ppf{end+1} = @(t) r*exp(1i*t) + mobius(0,a,R);
if v, Z = ppf{end}(t);  fill(real(Z),imag(Z),c); end

% external boundary
r = max(abs(C(1,:)+1i*C(2,:)))*1.05;
ppf{end+1} = @(t) r*exp(1i*t);
if v
    Z = ppf{end}(t); fillout(real(Z),imag(Z),[-1,1,-1,1]*r*1.1,c); 
    plot(real(Z),imag(Z),'k')
    hold off, axis equal, axis([-3,3,-3,3])
    axis off; %xticklabels(''); yticklabels('')
end

