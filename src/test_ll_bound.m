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

gam1 = repmat(u0,[1 J]); gam2 = repmat(sigma20, [1 J]);
alpha = rand([K J]); beta = rand([K J]);

%% Solve for log-likelihood bound

% Compute the log-likelihood bound
[ ll ] = ll_bound( r, n, u0, sigma20, M0, gam1, gam2, alpha, beta );


%% Plot the LL bound varying M0
v = logspace(-3 , 2, 50);
ll_M0 = arrayfun(@(x)ll_bound( r, n, u0, sigma20, x, gam1, gam2, alpha, beta ), v);
semilogx(v,ll_M0) 

%% Plot the LL bound varying u0
v = linspace(0 , 1, 50);
ll_u0 = arrayfun(@(x)ll_bound( r, n, x, sigma20, M0, gam1, gam2, alpha, beta ), v);
plot(v,ll_u0) 

%% Plot the LL bound varying sigma20
v = logspace(-3 , 2, 50);
ll_sigma20 = arrayfun(@(x)ll_bound( r, n, u0, x, M0, gam1, gam2, alpha, beta ), v);
semilogx(v,ll_sigma20) 

for i = 1:10
%% Open pool
% matlabpool open

%% Solve for u0
u0 = mean(gam1);

%% Solve for sigma20
sigma20 = mean(gam2 + (gam1-u0).^2);

%% Solve for M0
options = optimset('Display','iter');
% M0 = fminbnd(@(x)-ll_bound(r, n, u0, sigma20, x, gam1, gam2, alpha, beta),0,1e10)
M0 = fzero(@(x)dLdM0( x, u0, alpha, beta, gam1, gam2 ), M0, options);


%% Solve for alpha & beta
options = optimset('Display','off');
for j = 1:J
    for k = 1:K
         alpha(k,j) = fzero(@(x)dLdalpha(x, beta(k,j),M0,u0,r(k,j),n(k,j)), alpha(k,j), options);
    end
end

options = optimset('Display','off');
for j = 1:J
    for k = 1:K
         beta(k,j) = fzero(@(x)dLdbeta(alpha(k,j), x, M0,u0,r(k,j),n(k,j)), beta(k,j), options);
    end
end

%% Solve for gamma1

% gam1 = u0 - (sigma20*M0).*mean(psi(alpha)-psi(beta));
options = optimset('Display','off');
for j = 1:J
    gam1(j)= fzero(@(x)dLdgam1( x, gam2(j), alpha(:,j), beta(:,j), M0, u0, sigma20), gam1(j),...
        options);
end
plot(gam1)

%% Solve for gamma2
% options = optimset('Display','iter');
options = optimset('Display','off');
for j = 1:J
    gam2(j)= fzero(@(x)dLdgam2( gam1(j), x, M0, sigma20, J ), gam2(j),options);
end
plot(gam2)
disp('Done.')

end

%%

post_mu = gam1;
plot(theta(:,1:10)','.','MarkerSize',14)
hold on
plot(post_mu(1:10))
post_theta = alpha./(alpha+beta);
plot(post_theta(:,1:10)','o','MarkerSize',14)
hold off

%%
matlabpool close