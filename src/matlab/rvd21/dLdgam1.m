function [ f ] = dLdgam1( gam1, gam2, alpha, beta, M0, u0, sigma20)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[K,J] = size(alpha);

Eq1 = integral(@Eq1kernel,0,1);
Eq2 = integral(@Eq2kernel,0,1);

f = -(K/sigma20)*(gam1-u0) +M0*sum(psi(alpha)-psi(beta),1) -K*Eq1 -K*Eq2;

        
    function [f] = Eq1kernel(x)
        f = (1/gam2)*normpdf(x,gam1,sqrt(gam2)).*(x-gam1).*gammaln(x*M0);
    end

    function [f] = Eq2kernel(x)
        f = (1/gam2)*normpdf(x,gam1,sqrt(gam2)).*(x-gam1).*gammaln((1-x)*M0);
    end

end

