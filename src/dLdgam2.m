function [ f ] = dLdgam2( gam1, gam2, M0, sigma20, J )
%DLDgam2 Derivative of the bound with respect to gam2
%   Detailed explanation goes here

if gam2 <0
    gam2 = eps;
end

Eq1 = integral(@(x)Eq1kernel(x, M0, gam1, gam2), 0, 1);
Eq2 = integral(@(x)Eq2kernel(x, M0, gam1, gam2), 0, 1);

f = -J/(2*sigma20) - Eq1 -Eq2 + 1/(2*gam2);

    function [f] = Eq1kernel(x, M0, gam1, gam2)
        f = (2*gam2).^(-2).*normpdf(gam1,sqrt(gam2)).*(x-gam1).^2.*gammaln(x*M0) ...
                -pi*(2*pi*gam2)^(-3/2)*normpdf(gam1,sqrt(gam2)).*gammaln(x*M0);
    end
    function [f] = Eq2kernel(x, M0, gam1, gam2)
        f = (2*gam2).^(-2).*normpdf(gam1,sqrt(gam2)).*(x-gam1).^2.*gammaln((1-x)*M0) ...
                -pi*(2*pi*gam2)^(-3/2)*normpdf(gam1,sqrt(gam2)).*gammaln((1-x)*M0);
    end
end

