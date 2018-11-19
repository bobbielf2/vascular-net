function t_out = adaptive_panel2(f,p,t_in,tol,corner,stdpan)
% ADAPTIVE_PANEL2: panel choice given curve parameterization function(s)
%
% t_out = adaptive_panel2(f,p,t_in,tol,corner,stdpan,lam)
%
% Input:
%   f       Geometry parameterization func, or cell array of funcs, such as
%           Z: [a,b] -> C, speed: [a,b] -> R, giving curve in complex plane.
%   p       Use p-point Chebyshev nodes on each panel
%   t_in    Initial panel t_in, as 2-entry vector [a,b], parameter domain
%   tol     Prescribed tolerance
%   corner  Struct containing corner info
%           corner.tc = sorted list of corner locations in [a,b]
%           corner.lam = corner refinement rates corresp. to each corner.tc
%   stdpan  struct for standard panel, must have:
%           stdpan.nodes      - p nodes on [-1,1], col vec.
%           stdpan.testpts    - n test points on [-1,1], col vec.
%           stdpan.interpmat  - n*p interpolation matrix from values at nodes to
%                               values at test pts.
% Output:
%   t_out   row vec of parameter endpoints of the refined panels
%
% Notes:
% * removed globals (can't have these in any decent implementation!!), replaced
%   by new input needed - filled by below helper func
% * needs: gauss.m, interpmat_1d.m, adaptive_panel2_setup.m
% * m renamed p

% v2 by Alex Barnett 10/29/18; rewrite of adaptive_panel by Wu & Hai.

if nargin == 0, test_adaptive_panel; return, end
if isempty(corner), corner.tc = []; corner.lam = []; end
ncorn = numel(corner.tc);
nf = numel(f); fun = f; if nf==1, fun = {f}; end     % how many func handles?
n = numel(stdpan.testpts);
a = t_in(1); b = t_in(2); L = b-a;    % parameter length of this segment
i0 = find(corner.tc>a-L/2 & corner.tc<=a+L/2,1);  % corner near lower end?
i1 = find(corner.tc>a+L/2 & corner.tc<=b+L/2,1);  % corner near upper end?
split = [];    % empty will mean don't split
if L>tol
  if ~isempty(i0)      % corner-touching based splitting...
    split = a + L/corner.lam(i0);
  elseif ~isempty(i1)
    split = b - L/corner.lam(i1);
  else                 % or interp-based test for accuracte rep on panel...
    nodevals = nan(p,nf); testvals = nan(n,nf);
    for i=1:nf, nodevals(:,i) = fun{i}(a+L/2*(1+stdpan.nodes)); end
    for i=1:nf, testvals(:,i) = fun{i}(a+L/2*(1+stdpan.testpts)); end
    relerrs = (testvals - stdpan.interpmat*nodevals)./sum(abs(testvals));
    % TEST: approx arc length of panel, assume fun{1} is the curve parameterization
    arclen = fun{1}(a+L/2*(1+stdpan.nodes));
    arclen = sum(abs(arclen(2:end) - arclen(1:end-1)));
    if norm(relerrs,inf)>tol ...
            || max(L,arclen) > 2*pi/2/log10(1/tol) %bw - temporary rule restricting panel length (need split criterion for close target)
        split = (a+b)/2;       % symmetric split
    end
  end
end
if isempty(split)
  t_out = t_in;
else
  t0 = adaptive_panel2(f,p,[a split],tol,corner,stdpan);
  t1 = adaptive_panel2(f,p,[split b],tol,corner,stdpan);
  t_out = [t0(1:end-1) t1];      % remove duplicate endpt
end


%%%%%%%%%%%
function test_adaptive_panel

for tol = 10.^-(3:3:12)  % ========== tol loop
  figure;
  [p, stdpan] = adaptive_panel2_setup(tol);
  f = @(t) sin(1.5*t+1);    % test interpmat with around the right wiggliness
  fprintf('p=%d: interpmat err = %.3g\n', p, norm(stdpan.interpmat * f(stdpan.nodes) - f(stdpan.testpts),inf))
  xtest = -0.2 +0.3i;           % far interior pt in both domains
  
  for dom=1:2      % --------- loop over domains
    if dom==1  % smooth starfish.  set up Z, Z', Z''
      a = 0.3; w = 3; b = 0.5; Z = @(t) (1+a*cos(w*t + b)).*exp(1i*t);
      Zp = @(t) -w*a*sin(w*t + b).*exp(1i*t) + 1i*Z(t);
      Zpp = @(t) 1i*Zp(t) -w^2*a*cos(w*t + b).*exp(1i*t) -1i*w*a*sin(w*t + b).*exp(1i*t);
      corner = [];
    elseif dom==2  % Bremer boomerang on [-pi,pi]
      b = 3*pi/2; a = 2/tan(b/2)/3;                 % reentrant angle
      Z = @(t) -a*sin(3/2*abs(t))+1i*sin(t);        % t=0 is corner, for rel acc
      Zp = @(t) -a*3/2*sign(t).*cos(3/2*t)+1i*cos(t);
      Zpp = @(t) -a*(3/2)^2.*sin(3/2*abs(t))-1i*sin(t);
      corner.tc = 0; corner.lam = 3.0;   % for now
    end
    paramdomain = [-pi,pi];        % defines full curve
    funset = {Z, Zpp, @(t) abs(Zp(t))};   % for adapt, better than single func
    tpan = adaptive_panel2(funset,p,paramdomain,tol,corner,stdpan);   % do it
    npan = numel(tpan)-1; N=p*npan; s.x = nan(N,1); s.w=s.x; s.xp=s.x;
    % now build panel curve quadr...
    for m=1:npan       % for each panel, get param half-length (hl), nodes...
      hl=(tpan(m+1)-tpan(m))/2; t=tpan(m)+hl*(1+stdpan.nodes);
      i = (1:p)+(m-1)*p;   % indices in output arrays
      s.x(i) = Z(t); s.xp(i) = Zp(t);
      s.w(i) = hl*stdpan.weights';
    end
    cauchy = (1/(2i*pi))*sum(s.w.*s.xp./(xtest-s.x));  % test interior targ
    exact = -1.0; cauchyerr = abs(cauchy-exact);       % (cmpx density = 1)
    fprintf('dom=%d\t tol=%.3g\t p=%d\t npan=%d\t N=%d\terr=%.3g\n',dom,tol,p,npan,N,cauchyerr)
    subplot(1,2,dom);
    tt = linspace(paramdomain(1),paramdomain(2),1e3); plot(Z(tt),'-');
    hold on; plot(Z(tpan),'k+'); axis equal tight
    plot(s.x,'b.','markersize',8); plot(xtest,'rx');
    title(sprintf('\\epsilon=%.3g: p=%d, N=%d',tol,p,N))
  end              % -----------
  
end                     %  =========
