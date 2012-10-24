function [ f ] = dLdbeta( alpha, beta, M0, u0, r, n )
%DLDBETA Derivative of the bound with respect to beta
%   Detailed explanation goes here

if beta <= 0, beta = eps; end

f = psi(1,beta)*((1-u0)*M0-1 +(n-r) -beta) -psi(1,alpha+beta)*(M0+n-(alpha+beta));
end