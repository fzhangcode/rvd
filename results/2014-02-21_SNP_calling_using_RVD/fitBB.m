function [] = fitBB(runtimeResultsPath)
% Estimate the Beta-Binomial parameters using clinical data one lane at a
% time

% Author: Pat Flaherty
% Created:
% Modified: May 31, 2011
%% Fit the model to the reference data
load(fullfile(runtimeResultsPath, 'sampleInfo'),...
  'SampleInfo');
load(fullfile(runtimeResultsPath, 'seqAnnot'),...
  'SeqAnnot');
load(fullfile(runtimeResultsPath, 'filSeqData'),...
  'r','n');

nSample = length(SampleInfo);
isRef = SampleInfo.isReference;
nPos = length(SeqAnnot);

% Fit the Beta-Binomial model using the EM algorithm to the reference data
disp('Processing reference.')
[ M0_hat, mu0_hat, theta0_hat ] ...
	= beta_bino_em( r(:,isRef), n(:,isRef) );

% Store the model parameters in a structure
BBParamMLE0.M = M0_hat;
BBParamMLE0.mu = mu0_hat;
BBParamMLE0.theta = theta0_hat;

% Save the data 
save(fullfile(runtimeResultsPath,'BBParam'),'BBParamMLE0');