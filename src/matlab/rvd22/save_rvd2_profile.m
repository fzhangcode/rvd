function save_rvd2_profile(csvname)
% Make a read depth chart profile to test algorithm

rng(1000,'twister');

K = 2; J = 500;

u0 = 0.1; sigma20 = 0.01^2; M0 = 5e3; 

mu = normrnd(u0,sqrt(sigma20),[1 J]);
mu(mu<0) = 0; mu(mu>1) = 1;

theta = NaN(K,J);
for k = 1:K
    theta(k,:) = betarnd(mu.*M0,M0.*(1-mu)); 
end

n = repmat(1e5,K,J);
r = binornd(n,theta);

M = zeros(2*K,J);
M(1:2:2*K,:) = r;
M(2:2:2*K,:) = n;


csvwrite(csvname,M.')