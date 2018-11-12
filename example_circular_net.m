function [C, pp, ppf] = example_circular_net(smoothing)
%% generate a example (hydrodynamical) channel in a circular domain based on
% the topology of a "Fibonacci tree" (a fractal tree).
% Input:
%   smoothing   toggle: 
%                   0 = no smoothing, 
%                   1 = smoothed periodic cubic spline 
%                   2 = smoothed via partition of unity (test)
H = [6,7,6,5,5]; % num of layers for each main branch
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

addpath(genpath('fibchan'))
addpath(genpath('util'))

if nargin == 0, smoothing = 1; end
T = FibTree(H(end));
m = 5;
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
figure
c = [1,1,1]/2; % color
ind = find(isnan(C(1,:)));
pp = cell(numel(ind)-1,1);
ppf = cell(numel(ind)-1,1);
for i = 1:numel(ind)-1
    if smoothing == 1
        z = C(:,ind(i)+1:ind(i+1)-1); z = z(1,:)+1i*z(2,:);
        z = interp1(z, 1:.5:length(z)); % include midpoints of polygon edges, make interpolant closer to polygon.
        n = length(z);
        pp{i} = spcsp(1:n,z,0.98);
        t = linspace(0,2*pi,1000);
        ppf{i} = @(t) ppfun(pp{i},t);
        zz = ppf{i}(t);
        x = real(zz);
        y = imag(zz);
    elseif smoothing == 2
        z = C(:,ind(i)+1:ind(i+1)-2); z = z(1,:)+1i*z(2,:);
        t = linspace(0,2*pi,1000);
        ppf{i} = @(t) smooth_poly(z,t);
        zz = ppf{i}(t);
        x = real(zz);
        y = imag(zz);
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
        t = linspace(0,2*pi,1000);
        zz = ppf{i}(t);
        x = real(zz);
        y = imag(zz);
    end
    fill(x,y,c)
    hold on
end

% internal (source) boundary
r = 0.12;
t = linspace(0,2*pi,200);
ppf{end+1} = @(t) r*exp(1i*t) + mobius(0,a,R);
Z = ppf{end}(t);
fill(real(Z),imag(Z),c)

% external boundary
r = max(abs(C(1,:)+1i*C(2,:)))*1.03;
ppf{end+1} = @(t) r*exp(1i*t);
Z = ppf{end}(t);
fillout(real(Z),imag(Z),[-1,1,-1,1]*r*1.1,c);


hold off
axis equal
axis([-3,3,-3,3])
