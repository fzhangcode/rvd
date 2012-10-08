rng(1000,'twister');

K = 3; J = 1000;

u0 = 0.5; sigma20 = 0.1^2;

mu = normrnd(u0,sqrt(sigma20),[1 J]);

M0 = 1000; theta = NaN(K,J);
for k = 1:K
    theta(k,:) = betarnd(mu.*M0,M0.*(1-mu)); 
end

n = repmat(1000,K,J);
r = binornd(n,theta);

gam1 = repmat(0.5, [1 J]); gam2 = repmat(0.5, [1 J]);
alpha = repmat(3,[K J]); beta = repmat(2, [K J]);

% Compute the log-likelihood bound
[ ll ] = ll_bound( r, n, u0, sigma20, M0, gam1, gam2, alpha, beta );

% Solve for u0
u0 = mean(gam1) + (M0./(K*J)).*sum(sum(psi(alpha) - psi(beta)));

% Solve for sigma20

% Solve for M0
[ f ] = dLdM0( M0, u0, alpha, beta, gam1, gam2 );
M0 = fsolve( @(x)dLdM0(x, u0, alpha, beta, gam1, gam2), M0, ...
    'Display','iter' );

% Plot the LL bound varying M0
v = logspace(-1 , 3, 10);
ll_M0 = arrayfun(@(x)ll_bound( r, n, u0, sigma20, x, gam1, gam2, alpha, beta ), v);
plot(v,ll_M0) 

% Plot the LL bound varying u0
v = linspace(0 , 1, 100);
ll_u0 = arrayfun(@(x)ll_bound( r, n, x, sigma20, M0, gam1, gam2, alpha, beta ), v);
plot(v,ll_u0) 

% Plot the LL bound varying sigma20
v = logspace(-2 , 1, 100);
ll_sigma20 = arrayfun(@(x)ll_bound( r, n, u0, x, M0, gam1, gam2, alpha, beta ), v);
semilogx(v,ll_sigma20) 


