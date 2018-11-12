function z = smooth_poly(polygon,t)
% Smooth and parameterize a polygon based on partition of unity algorithm
% Input
%   polygon   complex coordinates of polygon vertices
%   t         query points of the parameterization, t in [0,2pi]
% Output
%   z         z = f(t) parameterization of the smoothed polygon on [0,2pi]
%
% written by Hai Zhu, annotated by Bowei Wu

polygon = preprocess(polygon);
z = assembly(polygon,t);
end

function v = preprocess(node)
% add a few points on each side

alpha = 2/(1+sqrt(5));
rec = [alpha^4; alpha^3; alpha^2; 1-alpha^2; 1-alpha^3; 1-alpha^4];

N = numel(node);
v = [];
for i = 1:N
    a = node(i); b = node(mod(i,N)+1);
    v = [v; a*(1-rec)+b*rec];
end
end

function fx = assembly(polygon,t)
% assemble smooth function G{i} between a{i}.x and a{i+1}.x
% 
NN = numel(polygon);
for i = 1:NN
    a{i}.x = polygon(i);
end
for i = 1:NN
    a{i}.gl = polynomial_basis(a{mod(i-2,NN)+1}.x,a{i}.x,a{mod(i,NN)+1}.x);
end

% calculate G{i} for given x
xx = t*NN/(2*pi);
fx = NaN(size(xx));
for k = 1:NN
    idx = find(xx>=(k-1) & xx<=k);
    xi = xx(idx)-(k-1);
    fx(idx) = pou(xi).*a{k}.gl(xi)+pou(1-xi).*a{mod(k,NN)+1}.gl(xi-1);
end

end

function Eta = pou(t)
% partition of unity
%
delta = 1/3; a = 1-2*delta;
h = @(s) exp(2*exp(-1./s)./(s-1));
Eta=h((t-delta)/a)./(h((t-delta)/a)+h(1-(t-delta)/a));
Eta(t>=0&t<=delta)=1;
Eta(t>=1-delta&t<=1)=0;
end



%%--------------------------------------------------------------------
function gl = polynomial_basis(xl,xc,xr)
% polynomial interpolation basis function at each node on interval [-1,1]
% 
% if flag
    gl = @(t) (1/2-t/2).^2.*xl + 2*(1/2+t/2).*(1/2-t/2).*xc + (1/2+t/2).^2.*xr;
% else
%     ql = @(t) -t.*xl + (t+1).*xc;   
%     qr = @(t) (1-t).*xc + t.*xr;
% 
%     gl = @(t) (1-t)/2.*ql(t) + (t+1)/2.*qr(t);
% end
end
