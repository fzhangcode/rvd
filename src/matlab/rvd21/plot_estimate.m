function plot_estimate( r, n, alpha, beta, gam1, u0, sigma20 )
%PLOT_ESTIMATE Plot the current posterior estimates for RVD2 model
%   Detailed explanation goes here
    clf
    rn = [r./n];
    thetaHat = alpha./(alpha+beta);
    
    plot(rn(:,1:10)','.','MarkerSize',14);
    hold on
    plot(thetaHat(:,1:10)','s','MarkerSize',8)
    plot(gam1(1:10),'--k')
    line([1 10],[u0 u0])
    line([1 10],[u0 u0]+sqrt(sigma20))
    line([1 10],[u0 u0]-sqrt(sigma20))
    hold off
end

