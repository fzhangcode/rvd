function [] = testBBAll(runtimeResultsPath, sDiff, alpha)
% Test for high error rates compared to null model for all samples

% Author: Pat Flaherty
% Created:
% Modified: October 19, 2012
%% Test observed 
% alpha = 1e-6;
load(fullfile(runtimeResultsPath,'seqAnnot.mat'),...
  'SeqAnnot'); % feature annotation
load(fullfile(runtimeResultsPath,'sampleInfo.mat'),...
  'SampleInfo'); % sample annotation
nSample = length(SampleInfo);

load(fullfile(runtimeResultsPath,'filSeqData.mat'),...
  'r','n','filDepth'); % sequence data
n0 = mean(n(:,SampleInfo.isReference),2);

load(fullfile(runtimeResultsPath,'BBParam.mat'),...
  'BBParamMLE0'); % model parameters
BBParamMLE0.M = BBParamMLE0.M(:,1);
BBParamMLE0.mu = BBParamMLE0.mu(:,1);
BBParamMLE0.theta = [];


for i = 1:nSample
	[ h1(:,i), p_val1(:,i), z_stat1(:,i), thetaBinMLE(:,i)] ...
		= testbb( BBParamMLE0.M, BBParamMLE0.mu, n0, r(:,i), n(:,i), alpha);
  
	deltaError1(:,i) = thetaBinMLE(:,i) - BBParamMLE0.mu;
  
  [ secondBase1(:,i) , secondBaseRate1(:,i)] ...
    = findsecondbase( filDepth(i), SeqAnnot.AlignRefBase );
end

h1plus = h1 & deltaError1 > sDiff;
%% Save the hypothesis test results
%
% Contents: 
% H: null: sample error rate = reference error rate (Binomial Test) true if
% reject.
% Hplus: H is rejected and the Binomial MLE is 0.1% greater than the 
% reference estimate of mu in the model
% HConfirm: H is true on both replicate lanes
% HPlusConfirm: Hplus is true on both replicat elanes
% SecondBase: the second most frequent base
% Second Base Rate: read depth for second base / total read depth
% Pvalue: Binomial test p-value
% Zstat: Z statistic for binomial test
% DeltaError: the difference between the estimated reference error rate and
% the estimated sample error rate
BBTest.H = h1;
BBTest.Hplus = h1plus; 
BBTest.SecondBase = secondBase1;
BBTest.SecondBaseRate = secondBaseRate1;
BBTest.Pvalue = p_val1;
BBTest.Zstat = z_stat1;
BBTest.DeltaError = deltaError1;

save(fullfile(runtimeResultsPath, 'BBTest.mat'), 'BBTest');