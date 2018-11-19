function [x, w, D, v] = cheby(N)
% based on chebyshev extreme cheb.m from spectral methods in matlab
% Chebyshev nodes, weights, and spectral differentiation matrix on [-1,1]
% 08/31/2016 Hai
% modified 3/20/18 B. Wu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Mathematical reference]:
%   Kuan Xu, "The Chebyshev points of the first kind", Applied Numerical
%   Mathematics 102 (2016) pp 17-30.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% chebyshev nodes
% theta = pi*(2*(1:N)'-1)/(2*N);
% x = -cos(theta);
theta2 = pi*(N - 2*(1:N)'+1)/(2*N);
x = -sin(theta2);
% better way of computing the nodes. Ref: Xu (sec 2.1, p 19)

% chebyshev weights
l = floor(N/2)+1;
K = 0:N-l;   
v = [2*exp(1i*pi*K/N)./(1-4*K.^2)  zeros(1,l)];
w = real(ifft(v(1:N) + conj(v(N+1:-1:2))))';
% j = 1:floor(N/2);
% w = 2/N*(1 - 2*sum( cos( repmat(j,N,1).*repmat((2*(0:N-1)'+1)*pi/N,1,floor(N/2)))...
%     .* (1./repmat(4*j.^2-1,N,1)), 2));

% spectral differentiation matrix
if nargout >= 3
    X = repmat(x,1,N);
    dX = X-X';  % x_i-x_j
    a = prod(dX+eye(N),2);  % cardinal function coeff
    D = (a*(1./a)')./(dX+eye(N));   % off diagonal element
    D = D - diag(sum(D,2)); % stable computation based on interpolation of constant function 1 (derivative 0)
end

% compute (true form) barycentric weights. Ref: Xu (sec 2.7, p 22)
if nargout >= 4
    v = ones(N, 1);
    v(end-1:-2:1) = -1;
    v = v.*cos(theta2);
end