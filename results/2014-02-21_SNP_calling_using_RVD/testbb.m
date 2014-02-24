function [ h, p_val, z_stat, varargout ] = testbb( M0_hat, mu0_hat, n0, r, n, alpha)
%TESTBB Normal test for Beta-Binomial Model
%   Tests for a significant difference between error rates under the
%   Beta-binomial model using a Gaussian approximation.

% Check that mu0_hat, n0, r, n are all same length
L = length(mu0_hat);
if(any([length(n0),length(r),length(n)] ~=L))
	error('Lengths of mu0_hat, n0, r, n must be equal.');
end

% Compute null distribution statistics
%M0_hat = M_hat(1); % prior read counts
%n0 = 50000; % reference read count
%mu0_hat = mu_hat(:,1); % reference error rate

% Null distribution standard deviation from Beta-Binomial model
sig0	= sqrt((mu0_hat.*(1-mu0_hat)./n0).*(1+(n0-1)./(M0_hat+1)));	

% Observed sample binomial error rate
theta_BinMLE = r./n; 

% Z-test for significant difference between r/n and reference r0/n0.
p_val	= 1-normcdf(theta_BinMLE, mu0_hat, sig0);
z_stat	= (theta_BinMLE - mu0_hat)./sig0;
h = p_val < alpha;

varargout{1} = theta_BinMLE;
varargout{2} = sig0;
end

