function plot_estimate( r, n, alpha, beta, gam1, gam2, u0, sigma20 )
%PLOT_ESTIMATE Plot the current posterior estimates for RVD2 model
%   Detailed explanation goes here
    N = 100;

    if length(gam2) == 1
        gam2 = repmat(gam2,size(gam1));
    end

    clf
    rn = [r./n];
    thetaHat = alpha./(alpha+beta);
    
    plot(rn(:,1:N)','.','MarkerSize',14);
    hold on
    plot(thetaHat(:,1:N)','s','MarkerSize',8)
    errorbar(gam1(1:N),sqrt(gam2(1:N)),'--k')
    line([1 N],[u0 u0])
    line([1 N],[u0 u0]+sqrt(sigma20))
    line([1 N],[u0 u0]-sqrt(sigma20))
    hold off
end

