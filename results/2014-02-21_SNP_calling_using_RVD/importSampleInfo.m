% Import and save sample information for clinical dataset

% Author: Pat Flaherty
% Created:
% Modified: November 7, 2012
function importSampleInfo(runtimeResultsPath, indexFileName)

% Import the sample info from the text file
fid = fopen(fullfile(indexFileName));
assert(fid >=3,'Could not open sample index file.');

C = textscan(fid, '%s%s%c%s%s', 'Delimiter',',', 'HeaderLines', 1);
fclose(fid);
[tag, sampleName, isRef, pair1BamFile, pair2BamFile] = C{:};

% Convert data types for sample info
isRef = logical(isRef=='Y');

% Define the number of samples
nSample = length(tag);

% Construct the dataset table
SampleInfo = dataset({tag,'Tag'},{sampleName,'SampleName'},...
{isRef,'isReference'},{pair1BamFile,'Pair1BamFile'},{pair2BamFile,'Pair2BamFile'});
SampleInfo = set(SampleInfo,'Description', 'Sample Information for Clinical Gene Data');

% Sort the sample info data
SampleInfo = sortrows(SampleInfo,{'SampleName'});

% Sort the sample table
save(fullfile(runtimeResultsPath,'sampleInfo.mat'),'SampleInfo','nSample');
end