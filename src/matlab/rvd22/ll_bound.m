function [ ll ] = ll_bound( r, n, u0, sigma20, M0, gam1, gam2, alpha, beta )
%LL_BOUND Variational bound on data log-likelihood
[K,J] = size(n);

EqlogPmu = -0.5*log(2*pi*sigma20) -(0.5/sigma20)*(gam2 + (gam1 - u0).^2);
EqlogPmu = repmat(EqlogPmu, [K 1]);

for j = 1:J
    EqlogG1(j) = integral(@(x)logG1(x,gam1(j),gam2,M0),0,1);
end
for j = 1:J
    EqlogG2(j) = integral(@(x)logG2(x,gam1(j),gam2,M0),0,1);
end

EqlogPtheta = gammaln(M0) -repmat(EqlogG1+EqlogG2,[K 1]) ...
    +(u0*M0-1)*(psi(alpha) - psi(alpha+beta)) ...
    +((1-u0)*M0-1)*(psi(beta) - psi(alpha+beta));

EqlogPr = gammaln(n+1) -gammaln(r+1) -gammaln(n-r+1) ...
    +r.*(psi(alpha)-psi(alpha+beta)) + (n-r).*(psi(beta)-psi(alpha+beta));

EqlogQmu = repmat(-0.5*log(2*pi*exp(1)*gam2),[1, J]);
EqlogQmu = repmat(EqlogQmu, [K 1]);

EqlogQtheta = gammaln(alpha+beta) - gammaln(alpha) - gammaln(beta) ...
    +(alpha-1).*psi(alpha) +(beta-1).*psi(beta) ...
    -(alpha+beta-2).*psi(alpha+beta);

ll = EqlogPmu + EqlogPtheta + EqlogPr - EqlogQmu - EqlogQtheta;
ll = sum(sum(ll));


end

function [f] = logG1(x, gam1, gam2, M0)
  f = normpdf(x,gam1,sqrt(gam2)).*gammaln(x*M0);
end
function [f] = logG2(x, gam1, gam2, M0)
  f = normpdf(x,gam1,sqrt(gam2)).*gammaln((1-x)*M0);
end