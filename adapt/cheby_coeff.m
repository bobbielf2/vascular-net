%% function CHEBY_COEFF to compute chebyshev interpolation coefficients
function [c,T] = cheby_coeff(fx)
% for the purpose of computing chebyshev interpolation coefficients
% get chebyshev polynomial coefficients
N = numel(fx);
c = zeros(size(fx));

% form coefficient computation matrix by Clenshaw (Numerical Methods for Special Functions)
T = 2*cos(pi*((1:N)-1)'*(2*(N:-1:1)-1)/(2*N))/N;
c(:) = T*fx(:);
end
