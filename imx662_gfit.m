function [p, perr] = imx662_gfit(x, y, p0, withBkg)
%IMX662_GFIT  Gaussian (+ optional constant background) fit with error bars.
%   p = [amplitude center sigma (background)].  perr = 1-sigma uncertainties
%   from the curvature (Hessian) of the chi-square at the minimum.
%   Uses base MATLAB only (fminsearch), Poisson-style weights.
    if nargin < 4, withBkg = false; end
    x = x(:); y = y(:);
    if withBkg
        f = @(p, x) p(1)*exp(-0.5*((x-p(2))./p(3)).^2) + p(4);
        if numel(p0) < 4, p0(4) = max(min(y), 0); end
    else
        f = @(p, x) p(1)*exp(-0.5*((x-p(2))./p(3)).^2);
        p0 = p0(1:3);
    end
    w = 1 ./ max(y, 1);
    chi2 = @(p) sum(w .* (y - f(p, x)).^2);
    opts = optimset('MaxFunEvals', 4e4, 'MaxIter', 4e4, 'TolX', 1e-9, 'TolFun', 1e-9, 'Display', 'off');
    p = fminsearch(chi2, p0, opts);
    p = fminsearch(chi2, p, opts);
    p(3) = abs(p(3));
    % numerical Hessian -> covariance = 2 * inv(H) * (reduced chi2)
    n = numel(p); H = zeros(n);
    h = max(abs(p) * 1e-4, 1e-6);
    for i = 1:n
        for j = i:n
            ei = zeros(1,n); ej = ei; ei(i) = h(i); ej(j) = h(j);
            H(i,j) = (chi2(p+ei+ej) - chi2(p+ei-ej) - chi2(p-ei+ej) + chi2(p-ei-ej)) / (4*h(i)*h(j));
            H(j,i) = H(i,j);
        end
    end
    dof = max(numel(x) - n, 1);
    red = chi2(p) / dof;
    C = 2 * pinv(H) * max(red, 1);
    perr = sqrt(abs(diag(C)))';
end
