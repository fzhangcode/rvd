function [ f ] = dLdgam2( gam1, gam2, M0, sigma20, J )
%DLDgam2 Derivative of the bound with respect to gam2
%   Detailed explanation goes here

if gam2 <0
    gam2 = eps;
end

Eq1 = integral(@(x)Eq1kernel(x, M0, gam1, gam2), 0, 1);
Eq2 = integral(@(x)Eq2kernel(x, M0, gam1, gam2), 0, 1);

f = -1/(2*sigma20) -(0.5/gam2)*Eq1 -(0.5/gam2)*Eq2 + 1/(2*gam2);

    function [f] = Eq1kernel(x, M0, gam1, gam2)
        f = (1/gam2).*normpdf(x,gam1, sqrt(gam2)).*gammaln(x*M0).*(x-gam1).^2 ...
            -normpdf(x, gam1,sqrt(gam2)).*gammaln(x.*M0);
            
    end
    function [f] = Eq2kernel(x, M0, gam1, gam2)
        f = (1/gam2).*normpdf(x, gam1, sqrt(gam2)).*gammaln((1-x)*M0).*(x-gam1).^2 ...
            -normpdf(x,gam1,sqrt(gam2)).*gammaln((1-x).*M0);    
    end
end

