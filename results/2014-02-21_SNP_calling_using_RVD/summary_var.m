function summary_var(resultsPath,runtimeResultsPath)
% Summarize the distribution statistics for the clinical data

% Author: Pat Flaherty
% Created:
% Modified: October,24 2012
%% Load the depth chart data
load(fullfile(runtimeResultsPath,'filSeqData.mat'));
load(fullfile(runtimeResultsPath,'seqAnnot.mat'));
load(fullfile(runtimeResultsPath,'sampleInfo.mat'));

%% Compute position-conditional and position-marginal variance for lane 1
fname = fullfile(resultsPath,'summary_error.txt');
fid = fopen(fname, 'wt');
fprintf(fid,'\n');

fprintf(fid,'Lane Summary Statistics\n');
isRef = SampleInfo.isReference;
theta = r(:,isRef)./n(:,isRef);

Vartheta = var(theta(:),1);
fprintf(fid,'Total Variation = %0.4e\n', Vartheta);

VarEtheta = var(mean(r(:,isRef)./n(:,isRef),2),1);
fprintf(fid,'Position-marginal variation (var): Var[E(r/n|pos)] = %0.4e\n',...
	VarEtheta);

EVartheta = mean(var(r(:,isRef)'./n(:,isRef)',1));
fprintf(fid,'Position-conditional variation (var): E[Var(r/n|pos)] = %0.4e\n',...
	EVartheta);

F = VarEtheta/EVartheta;

K = size(theta,1); % number of groups
N = numel(theta); % total number of samples

pval1 = 1-fcdf(F,K,N-K);

fprintf(fid,'H0: Equality between marginal and conditional p = %0.3e\n', pval1);

%% Display variance calculations
fprintf(fid,'\n');
fprintf(fid,'Average error rate: %0.4f\n', ...
  mean(sum(r(:,isRef))./sum(n(:,isRef))) );

fprintf(fid,'Average standard deviation: %0.3e\n',...
  sqrt(mean(var(r(:,isRef)./n(:,isRef)))) );

fprintf(fid,'Percentiles (5, 95): (%0.4f, %0.4f)\n',...
  mean(prctile(r(:,isRef)./n(:,isRef),[5 95]),2) );

fprintf(fid,'Min and Max Error rates: (%0.4f, %0.4f)\n',...
  [min(min(r(:,isRef)./n(:,isRef))), max(max(r(:,isRef)./n(:,isRef)))]);

fprintf(fid,'\n');
fclose(fid);