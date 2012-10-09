rng(1000,'twister');

K = 3; J = 100;

u0 = 0.5; sigma20 = 0.1^2;

mu = normrnd(u0,sqrt(sigma20),[1 J]);

M0 = 10; theta = NaN(K,J);
for k = 1:K
    theta(k,:) = betarnd(mu.*M0,M0.*(1-mu)); 
end

n = repmat(100,K,J);
r = binornd(n,theta);

gam1 = repmat(u0, [1 J]); gam2 = repmat(sigma20, [1 J]);
alpha = repmat(3,[K J]); beta = repmat(2, [K J]);

% Compute the log-likelihood bound
[ ll ] = ll_bound( r, n, u0, sigma20, M0, gam1, gam2, alpha, beta );

% Plot the LL bound varying M0
v = logspace(-3 , 2, 50);
ll_M0 = arrayfun(@(x)ll_bound( r, n, u0, sigma20, x, gam1, gam2, alpha, beta ), v);
semilogx(v,ll_M0) 

% Plot the LL bound varying u0
v = linspace(0 , 1, 50);
ll_u0 = arrayfun(@(x)ll_bound( r, n, x, sigma20, M0, gam1, gam2, alpha, beta ), v);
plot(v,ll_u0) 

% Plot the LL bound varying sigma20
v = logspace(-3 , 0, 50);
ll_sigma20 = arrayfun(@(x)ll_bound( r, n, u0, x, M0, gam1, gam2, alpha, beta ), v);
semilogx(v,ll_sigma20) 


