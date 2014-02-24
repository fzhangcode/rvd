function RVD(configFileName) 
% RVD Process depth charts to variant calls using the RVD algorithm
%
% Syntax:  RVD(configFileName)
%
% Inputs:
%    configFileName - Full path to configuration .txt file
%
% Outputs:
%    none
%
% Example: 
%    none
%
% Other m-files required: testBBAll, testbb, summary_var, importSequenceData,
% importSampleInfo, getRVDConfig, fitBB, findsecondbase, filterStrandError,
% filter_strand_by_diff, exportCallTable, createDepthChart, compileRawDepth,
% combineRawReads, beta_bino_em, bb_loglik_comp 
%
% Subfunctions: none 
% MAT-files required: none
%
% See also: none

% Author: Patrick Flaherty
% Copyright (c) 2012 Stanford University. All rights reserved.
% email: 
% Website: http://www.
% November 2012; Last revision: 07-November-2012

% Setup the environment
[rootPath, dataPath, resultsPath, ...
    sampleFile, ROI, minimumError, referenceFile, ...
    runtimeResultsPath, phredCutoff, qualOffset] = getRVDConfig(configFileName);
  
% Import clinical data
importSampleInfo(runtimeResultsPath, sampleFile);

% Import sequence data
importSequenceData(dataPath,runtimeResultsPath, ROI, referenceFile, phredCutoff, qualOffset);

% Preprocess and summarized clinical data
filterStrandError(runtimeResultsPath, 2.5e-3, 1e-6);
summary_var(resultsPath,runtimeResultsPath)

fprintf('\nBeginning model estimation...\n');
% Fit beta-binomial model to data
fitBB(runtimeResultsPath);

% Display summary of data
testBBAll(runtimeResultsPath, 1e-3, 1e-6);
exportCallTable(runtimeResultsPath, resultsPath);

fprintf('\nDone.\n');
end