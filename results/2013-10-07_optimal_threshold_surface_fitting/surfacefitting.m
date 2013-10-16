cov = h5read('optimalT.hdf5','/cov');
T_ED = h5read('optimalT.hdf5','/T_ED');
T_L1 = h5read('optimalT.hdf5','/T_L1');
T_MCC = h5read('optimalT.hdf5','/T_MCC');
dilution = 0.01*h5read('optimalT.hdf5','/dilutionList');


% average 8 times sampling
T_ED=mean(T_ED,3);
T_L1=mean(T_L1,3);
T_MCC=mean(T_MCC,3);

% log surface fitting
cov_log=arrayfun(@(x) log10(x), cov);
dilution_log=arrayfun(@(x) log10(x), dilution);

[fitresult1, gof1] = createFit_linear(dilution, cov, T_MCC);
[fitresult2, gof2] = createFit_log(dilution_log, cov_log, T_MCC);


optimalT2=@(x,y) fitresult2.p00+fitresult2.p10*x+fitresult2.p01*y+...
    fitresult2.p20*x.^2+fitresult2.p11*x.*y+fitresult2.p02*y.^2;

[xdilution_log, xcov_log]=meshgrid(dilution_log,cov_log);

optimalT=optimalT2(xdilution_log,xcov_log);

if ~exist('optimalT_logfit.hdf5','file')
    h5create('optimalT_logfit.hdf5','/T_fit',size(optimalT));
    h5write('optimalT_logfit.hdf5','/T_fit',optimalT);
end
