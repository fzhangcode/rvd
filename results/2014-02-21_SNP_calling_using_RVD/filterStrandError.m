function filterStrandError(runtimeResultsPath, sDiff, alpha)
% Filter strands based on error rate differential.
% Input:  runtimeResultsPath
%         sDiff = differential error rate cutoff
%         alpha = differential error rate p-value cutoff

% Author: Pat Flaherty
% Created: November 21, 2010
% Modified: March 23, 2012
%% Load the depth chart data for the dilution series
load(fullfile(runtimeResultsPath, 'sampleInfo.mat'));
load(fullfile(runtimeResultsPath, 'rawSeqData.mat'),...
  'errDepth', 'refDepth', 'rawDepth');

isRef = SampleInfo.isReference;

%% Filer for strand errors
[ r, n, errDepth, refDepth, filDepth] ...
  = combineRawReads( errDepth, refDepth, rawDepth, isRef, sDiff, alpha );

%% Export the filtered data 
save(fullfile(runtimeResultsPath, 'filSeqData.mat'),...
  'r', 'n', 'errDepth', 'refDepth','filDepth');