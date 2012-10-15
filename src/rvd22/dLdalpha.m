function [ f ] = dLdalpha( alpha, beta, M0, u0, r, n )
%DLDALPHA Derivative of the bound with respect to alpha
%   Detailed explanation goes here

if alpha <=0
    alpha = eps;
end

f = psi(1,alpha)*(u0*M0 +r -alpha) -psi(1,alpha+beta)*(M0+n-(alpha+beta));
end