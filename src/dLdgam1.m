function [ f ] = dLdgam1( gam1, gam2, alpha, beta, M0, u0, sigma20)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[K,J] = size(alpha);

Eq1 = integral(@Eq1kernel,0,1);
Eq2 = integral(@Eq2kernel,0,1);
Eq3 = integral(@Eq3kernel,0,1);
Eq4 = integral(@Eq4kernel,0,1);

f = -(K/sigma20)*(gam1-u0) -(1/gam2)*Eq1 +(gam1/gam2)*Eq2 ...
    - (1/gam2)*Eq3 + (gam1/gam2)*Eq4 + M0*sum(psi(alpha)-psi(beta),1);

        
    function [f] = Eq1kernel(x)
        f = normpdf(x,gam1,sqrt(gam2)).*x.*gammaln(x*M0);
    end

    function [f] = Eq2kernel(x)
        f = normpdf(x,gam1,sqrt(gam2)).*gammaln(x*M0);
    end

    function [f] = Eq3kernel(x)
        f = normpdf(x,gam1,sqrt(gam2)).*x.*gammaln((1-x)*M0);
    end

    function [f] = Eq4kernel(x)
        f = normpdf(x,gam1,sqrt(gam2)).*gammaln((1-x)*M0);
    end
end

