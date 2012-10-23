function [ll] = plot_profile(r, n, u0, sigma20, M0, alpha, beta, gam1, gam2)
% PLOT_PROFILE Plot the log-likelihood bound as a function of key
% parameters

% Plot the LL bound varying M0
v = logspace(-3 , 2, 50);
ll_M0 = arrayfun(@(x)ll_bound( r, n, u0, sigma20, x, gam1, gam2, alpha, beta ), v);
subplot(3,1,1), semilogx(v,ll_M0)
title('M_0')

% Plot the LL bound varying u0
v = linspace(0 , 1, 50);
ll_u0 = arrayfun(@(x)ll_bound( r, n, x, sigma20, M0, gam1, gam2, alpha, beta ), v);
subplot(3,1,2), plot(v,ll_u0)
title('u_0')

% Plot the LL bound varying sigma20
v = logspace(-3 , 2, 50);
ll_sigma20 = arrayfun(@(x)ll_bound( r, n, u0, x, M0, gam1, gam2, alpha, beta ), v);
subplot(3,1,3), semilogx(v,ll_sigma20)
title('\sigma_2^0')

% Compute the log-likelihood bound
[ ll ] = ll_bound( r, n, u0, sigma20, M0, gam1, gam2, alpha, beta );

end
