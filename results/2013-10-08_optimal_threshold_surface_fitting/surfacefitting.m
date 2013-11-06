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

%% abandon the 10% dilution rate
dilution3 = dilution(1:3);
T_ED3 = T_ED(:,1:3);
T_L13 = T_L1(:,1:3);
T_MCC3 = T_MCC(:,1:3);
% log surface fitting
cov_log3=cov_log;
dilution_log3=arrayfun(@(x) log10(x), dilution);


[fitresult311, gof] = poly311(dilution_log3, cov_log3, T_MCC);

optimalT311=@(x,y) fitresult311.p00 + fitresult311.p10*x + fitresult311.p01*y

[xdilution_log3, xcov_log3]=meshgrid(dilution_log3,cov_log3);

optimalT311=optimalT311(xdilution_log3,xcov_log3);

if ~exist('optimalT311.hdf5','file')
    h5create('optimalT311.hdf5','/T_fit',size(optimalT311));
    h5write('optimalT311.hdf5','/T_fit',optimalT311);
end
