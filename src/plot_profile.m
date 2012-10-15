%% Solve for log-likelihood bound

% Compute the log-likelihood bound
% [ ll ] = ll_bound( r, n, u0, sigma20, M0, gam1, gam2, alpha, beta );


%% Plot the LL bound varying M0
% v = logspace(-3 , 2, 50);
% ll_M0 = arrayfun(@(x)ll_bound( r, n, u0, sigma20, x, gam1, gam2, alpha, beta ), v);
% semilogx(v,ll_M0) 

%% Plot the LL bound varying u0
% v = linspace(0 , 1, 50);
% ll_u0 = arrayfun(@(x)ll_bound( r, n, x, sigma20, M0, gam1, gam2, alpha, beta ), v);
% plot(v,ll_u0) 

%% Plot the LL bound varying sigma20
% v = logspace(-3 , 2, 50);
% ll_sigma20 = arrayfun(@(x)ll_bound( r, n, u0, x, M0, gam1, gam2, alpha, beta ), v);
% semilogx(v,ll_sigma20) 

