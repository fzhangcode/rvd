function [rootPath, dataPath, resultsDir, ...
    sampleFile, ROI, minimumError, referenceFile, ...
    runtimeResultsPath, phredCutoff, qualOffset] = getRVDConfig(configFileName)
% Load configuration file for RVD
%
% Syntax:  getRVDConfig(configFileName)
%
% Inputs:
%    configFileName - Full path to configuration .txt file
%
% Outputs:
%    rootPath - current deployed path
%    dataPath - path to data directory
%    resultsDir - path to directory storing final call tables
%    sampleFile - full path to sample meta info csv file
%    ROI - region of interest for algorithm
%    minimumError - error resolution threshold for algorithm
%    referenceFile - fasta format reference sequence
%    runtimeResultsPath - directory to store temporary data
%    phredCutoff - base quality threshold for filtering low quality reads
%    ASCIIQualOffset (optional) - the ASCII offset used in determining
%        phred base quality. Default is 33.    
%
% Example: 
%    none
%
% Other m-files required: none 
% Subfunctions: none 
% MAT-files required: none
%
% See also: none

% Author: Patrick Flaherty
% Copyright (c) 2012 Stanford University. All rights reserved.
% email: 
% Website: http://www.
% November 2012; Last revision: 07-November-2012

% Define the project root path
% Edit this for your own installation.
if isdeployed
  rootPath = ctfroot;
  rootPath = fullfile(rootPath,'RVD');
else
  rootPath = pwd;
end


% Define the current results directory path
fid = fopen(fullfile(configFileName),'r');
assert(fid >=3, 'Could not open sample configuration file.');

dataPath      = fgetl(fid);
resultsDir    = fgetl(fid);
sampleFile    = fgetl(fid);
ROI           = fgetl(fid);
minimumError  = fgetl(fid);
referenceFile = fgetl(fid);
phredCutoff   = fgetl(fid);
qualOffset    = fgetl(fid);
if ischar(qualOffset)
    qualOffset = str2double(qualOffset);
else
    qualOffset = 33;
end
fclose(fid);

% Report configuration parameters to user to stdout
fprintf('\n');
fprintf('Data path is %s\n',dataPath);
fprintf('Results path is %s\n',resultsDir);
fprintf('Sample contents file is: %s\n',sampleFile);
fprintf('Region of Interest is %s\n',ROI);
fprintf('Allele frequency detection threshold is %s\n',minimumError);
fprintf('Reference sequence file is %s\n',referenceFile);
fprintf('Phred base quality score threshold is %s\n',phredCutoff);
fprintf('ASCII offset for base quality is %d\n',qualOffset);
% Create a results directory and define path vars
runtimeResultsPath = fullfile(resultsDir, 'algorithmRuntimeResults');
if ~exist(runtimeResultsPath,'dir')
  mkdir(runtimeResultsPath);
end

