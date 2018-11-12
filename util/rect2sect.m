function [xout,yout] = rect2sect(x,y,th,ifRnd,r,R)
% RECT2SECT - conformally map a rectangular region to a circular sector;
% can add random smooth distortion.
%
% Input:  x = tree height (normalized  0 < x < 1)
%         y = tree width (map to a circle without overlapping requires max(y) - min(y) < 2*pi)
%         th = optional angle rotated about origin
%         ifRnd = toggle for random distortion
%         r = radius of center circle
%         R = radius of network
%
% Requires: randfun.m
% Bowei Wu, 10-15-2018

if nargin < 6, R = exp(1); end
if nargin < 5, r = 0.25; end
if nargin < 4, ifRnd = 0; end
if nargin < 3, th = 0; end
% circular sector with random fluctuations
a = log(r); b = log(R);
xs = (b-a)*x.^.5+a;
if ifRnd
    % wave = @(x) 0.1.*sin(6*pi*x.^1.2);
    wave = @(x,y) 0.02*randfun(x.^1.1+0.05*nthroot(y/pi*3,3)+.1,5);
    z = exp(xs + 1i*(y + wave(x,y)) + 1i*th);
    %z = exp(xs + wave(x*.04,y*.04)*2 + 1i*(y + wave(x,y)) + 1i*th);
else
    z = exp(xs + 1i*y + 1i*th);
end
xout = real(z);
yout = imag(z);

end

