function [ f ] = dLdgam2( gam1, gam2, M0, sigma20 )
%DLDgam2 Derivative of the bound with respect to gam2
%   Detailed explanation goes here

J = length(gam1);

if gam2 <0
    gam2 = eps;
end

Eq1 = 0; Eq2 = 0;
for j = 1:J
    Eq1 = Eq1 + integral(@(x)Eq1kernel(x, M0, gam1(j), gam2), 0, 1);
    Eq2 = Eq2 + integral(@(x)Eq2kernel(x, M0, gam1(j), gam2), 0, 1);
end


f = -J/(2*sigma20) -(0.5/gam2)*Eq1 -(0.5/gam2)*Eq2 + J/(2*gam2);

    function [f] = Eq1kernel(x, M0, gam1, gam2)
        f = (1/gam2).*normpdf(x,gam1, sqrt(gam2)).*gammaln(x*M0).*(x-gam1).^2 ...
            -normpdf(x, gam1,sqrt(gam2)).*gammaln(x.*M0);
            
    end
    function [f] = Eq2kernel(x, M0, gam1, gam2)
        f = (1/gam2).*normpdf(x, gam1, sqrt(gam2)).*gammaln((1-x)*M0).*(x-gam1).^2 ...
            -normpdf(x,gam1,sqrt(gam2)).*gammaln((1-x).*M0);    
    end
end

