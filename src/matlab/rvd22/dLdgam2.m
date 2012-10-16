function [ f ] = dLdgam2( gam1, gam2, M0, sigma20 )
%DLDgam2 Derivative of the bound with respect to gam2
%   Detailed explanation goes here

J = length(gam1);

Eq1 = zeros(1,J); Eq2 = zeros(1,J);
for j = 1:J
    Eq1(j) = integral(@(x)Eq1kernel(x, M0, gam1(j), gam2), 0, 1);
    Eq2(j) = integral(@(x)Eq2kernel(x, M0, gam1(j), gam2), 0, 1);
end

f = -1/(2*sigma20) -(0.5/gam2)*mean(Eq1) -(0.5/gam2)*mean(Eq2) + 1/(2*gam2);

    function [f] = Eq1kernel(x, M0, gam1, gam2)
        f = (1/gam2).*normpdf(x,gam1, sqrt(gam2)).*gammaln(x*M0).*(x-gam1).^2 ...
            -normpdf(x, gam1,sqrt(gam2)).*gammaln(x.*M0);
            
    end
    function [f] = Eq2kernel(x, M0, gam1, gam2)
        f = (1/gam2).*normpdf(x, gam1, sqrt(gam2)).*gammaln((1-x)*M0).*(x-gam1).^2 ...
            -normpdf(x,gam1,sqrt(gam2)).*gammaln((1-x).*M0);    
    end
end

