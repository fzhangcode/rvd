%% Problem Setup

save_rvd2_profile('test_prof1.csv');
M = csvread('test_prof1.csv');
r = M(:, 1:2:size(M,2)).';
n = M(:, 2:2:size(M,2)).';
[K,J] = size(r);
disp('Generated test data.')
%% Estimate Model Parameters
u0 = 0.5; sigma20 = 0.01^2; M0=5000;

% plot_estimate( r, n, alpha, beta, gam1, u0, sigma20 )

% matlabpool open
[ u0,sigma20,M0,alpha,beta,gam1,gam2 ] ...
    = rvd2_est( r, n, u0, sigma20, M0 );
% matlabpool close

