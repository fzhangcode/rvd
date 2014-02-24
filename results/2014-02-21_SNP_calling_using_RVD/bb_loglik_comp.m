function [ ll] = bb_loglik_comp(mu, M, n, r, theta )
%BB_LOGLIK_COMP Beta-Binomial complete-data log-likelihood



[~,N] = size(r);
rep_mu = repmat(mu,[1 N]);

ll = gammaln(n+1) - gammaln(r+1) - gammaln(n-r+1) ...
	+ r.*log(theta) + (n-r).*log(1-theta) ...
	+ gammaln(M) - gammaln(rep_mu.*M) - gammaln((1-rep_mu).*M) ...
	+ (rep_mu.*M-1).*log(theta) + ((1-rep_mu).*M-1).*log(1-theta);

end

