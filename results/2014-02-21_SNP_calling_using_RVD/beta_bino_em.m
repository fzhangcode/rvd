function [ M_hat, mu_hat, theta_hat ] = beta_bino_em( r, n )
%BETA_BINO_EM Expectation-Maximization Algorithm for Beta-Binomial
%   EM algorithm for beta-binomial model. Inputs are the error read counts
%   (r) and the total number of reads (n). Outputs are the prior count
%   estimate (M), the position specific error rate estimate (mu) and the
%   position x replicate error rate (theta). r and n are matrices with
%   positions in rows and replicates in columns.
% Author: Pat Flaherty
% Created: February 1, 2011
% Modified: June 1, 2011

[J,N] = size(r); % read position x replicate

% Set the initial parameter values at their method-of-moments estimates
% 0 = wild-type, 1 = variant
mu_hat  = sum(r,2)./sum(n,2);	     	% prior mutant fraction
M_hat   = mean(mean(n));		% prior sequence depth

lc = -realmax; delta_lc=Inf;
lc_save = [];

% EM Algorithm Loop
% NOTE: From a code profiling, the majority of the execution time is in the
% fmincon functon. A substantial speed-up can be achieved by using the
% results and Blei and Minka.
while(delta_lc >1e-4)
	% E-step: estimate theta
	theta_hat = (repmat(mu_hat.*M_hat,[1 N]) + r -1)./(M_hat + n -2);
	
	% M-step: estimate {mu, M} using MLE
	
	% Optimization problem setup for M_hat
	problem_M.objective = @(x)(-sum(sum(bb_loglik_comp(mu_hat, x, n, r, theta_hat )) ));
	problem_M.x0 = M_hat;
	problem_M.lb = 0;
	problem_M.ub = Inf;
	problem_M.solver='fmincon';
	problem_M.options = optimset('Algorithm','interior-point',...
		'Display','off');
	
	[M_hat, ~, exitflag] = fmincon(problem_M);
	if(exitflag < 0), disp('Problem estimating M.'); end
	
	for j = 1:J
		% Optimization problem setup for mu_hat
		problem_mu.objective = @(x)(-sum(sum(bb_loglik_comp(x, M_hat, n(j,:), r(j,:),...
			theta_hat(j,:))) ));
		problem_mu.x0 = mu_hat(j);
		problem_mu.lb = 0;
		problem_mu.ub = 1;
		problem_mu.solver='fmincon';
		problem_mu.options = optimset('Algorithm','interior-point',...
			'Display','off');
		
		[mu_hat(j),~,exitflag] = fmincon(problem_mu);
		if(exitflag < 0), disp('Problem estimating mu.'); end
	end
	
	
	% Update the log-likehihood
	lc_old = lc;
	lc = sum(sum(bb_loglik_comp(mu_hat, M_hat, n, r, theta_hat )));
	delta_lc = (lc-lc_old)/abs(lc_old);
	lc_save = [lc_save lc];
	if(delta_lc <0), disp('Warning: Log-likelhood decreased.'); end
% 	plot(lc_save,'.-', 'LineWidth', 2), pause(0.5);
	
end

end

