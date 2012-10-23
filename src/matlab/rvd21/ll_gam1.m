function [ ll ] = ll_gam1( u0,sigma20,M0,gam1j,gam2j,alphaj,betaj )
%LL_GAM1 Log-likelihood bound with respect to j-th component of gam1
%   Detailed explanation goes here

K = length(alphaj);

EqlogG1j = integral(@(x)logG1(x,gam1j,gam2j,M0),0,1);
EqlogG2j = integral(@(x)logG2(x,gam1j,gam2j,M0),0,1);

ll = -0.5*(K/sigma20)*(gam1j - u0)^2 ...
        +gam1j*M0*sum(psi(alphaj)-psi(betaj)) ...
        -K*EqlogG1j -K*EqlogG2j;

    function [f] = logG1(x, gam1, gam2, M0)
        f = normpdf(x,gam1,sqrt(gam2)).*gammaln(x*M0);
    end
    function [f] = logG2(x, gam1, gam2, M0)
        f = normpdf(x,gam1,sqrt(gam2)).*gammaln((1-x)*M0);
    end

end

