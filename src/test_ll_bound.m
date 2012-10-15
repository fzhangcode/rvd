%% Problem Setup

rng(1000,'twister');

K = 3; J = 100;

u0 = 0.5; sigma20 = 0.1^2;

mu = normrnd(u0,sqrt(sigma20),[1 J]);
mu(mu<0) = 0; mu(mu>1) = 1;

M0 = 50; theta = NaN(K,J);
for k = 1:K
    theta(k,:) = betarnd(mu.*M0,M0.*(1-mu)); 
end

n = repmat(1000,K,J);
r = binornd(n,theta);

%% Estimate Model Parameters
matlabpool open
[ u0,sigma20,M0,alpha,beta,gam1,gam2 ] ...
    = rvd2_est( r, n, u0, sigma20, M0 );
matlabpool close

plot_profile(r,n,u0,sigma20,M0,alpha,beta,gam1,gam2)
