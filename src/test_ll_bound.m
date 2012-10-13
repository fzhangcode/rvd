%% Problem Setup

rng(1000,'twister');

K = 3; J = 100;

u0 = 0.5; sigma20 = 0.1^2;

mu = normrnd(u0,sqrt(sigma20),[1 J]);
sum(mu<0 |mu>1)
mu(mu<0) = 0; mu(mu>1) = 1;

M0 = 50; theta = NaN(K,J);
for k = 1:K
    theta(k,:) = betarnd(mu.*M0,M0.*(1-mu)); 
end

n = repmat(100,K,J);
r = binornd(n,theta);

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

%% Estimate Model Parameters
% matlabpool open
[ u0,sigma20,M0,alpha,beta,gam1,gam2 ] ...
    = rvd2_est( r, n, u0, sigma20, M0 );
% matlabpool close
%%

% post_mu = gam1;
% plot(theta(:,1:10)','.','MarkerSize',14)
% hold on
% plot(post_mu(1:10))
% post_theta = alpha./(alpha+beta);
% plot(post_theta(:,1:10)','o','MarkerSize',14)
% hold off
