function [] = exportCallTable(runtimeResultsPath, resultsPath)
% Exports call tables for BB model

% Author: Patrick Flaherty
% Created: May 31, 2011
% Modified: November 7, 2012

load(fullfile(runtimeResultsPath,'seqAnnot.mat'),...
  'SeqAnnot'); % feature annotation
load(fullfile(runtimeResultsPath,'sampleInfo.mat'),...
  'SampleInfo'); % sample annotation
nSample = length(SampleInfo);

load(fullfile(runtimeResultsPath,'filSeqData.mat'),...
  'filDepth'); % sequence data
load(fullfile(runtimeResultsPath,'BBParam.mat')); % model parameters
load(fullfile(runtimeResultsPath,'BBTest.mat')); % test data

%% Find Second Base and export to call table
 for i = 1:nSample
   Ha = BBTest.H(:,i);   
   [ secondBase , secondBaseRate] ...
     = findsecondbase( filDepth(i), SeqAnnot.AlignRefBase );
   [mu0] = BBParamMLE0.mu;

  % Write testing output to a table
	fname = sprintf('calltbl_%d_%3s.txt',i,SampleInfo.SampleName{i});
	fname = fullfile(resultsPath,fname);
  
  dset = dataset({SeqAnnot.AlignRefPos,'AlignReferencePosition'},...
    {SeqAnnot.AlignRefChar,'AlignmentBase'},...
    {Ha,'Call'}, ...
    {secondBase,'SecondBase'},...
    {BBTest.DeltaError(:,i).*100,'CenteredErrorPrc'},...
    {mu0.*100,'ReferenceErrorPrc'},...
    {secondBaseRate.*100,'SecondBasePrc'});
  
  export(dset,'file',fname);
end