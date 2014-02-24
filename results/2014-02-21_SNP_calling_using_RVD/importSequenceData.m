function importSequenceData(dataPath, runtimeResultsPath, ROI, referenceFile, cutoff_phred, qualOffset)
% Convert BAM file to data matrix for statistical analysis

% Author: Pat Flaherty
% Created:
% Modified: November 7, 2012

% Load the sample manifest data
load(fullfile(runtimeResultsPath, 'sampleInfo.mat'),...
  'SampleInfo','nSample')

%import the reference sequence
[Header, refSeq] = fastaread(referenceFile);

cutoff_phred = str2double(cutoff_phred);

% Import the depth chart data
rawDepth = struct([]);
fprintf('Conversting BAM file to depth charts...\n')
firstPos = 1;
lastPos = length(refSeq);
inBounds = zeros(4,nSample);
for iSample = 1:nSample
  fprintf('\tSample %d of %d. \n',iSample,nSample)
  % Create Pair 1 depth chart from BAM and read into p1Fwd and p1Rev
  f1Name = [SampleInfo.Pair1BamFile{iSample}];
  [p1Fwd, p1Rev, inBounds(1,iSample), firstPos, lastPos, inBounds(2,iSample)]...
      = createDepthChart(runtimeResultsPath, dataPath, refSeq, f1Name, cutoff_phred, firstPos, lastPos, qualOffset);

  % Create Pair 2 depth chart from BAM and read into p2Fwd and p2Rev
  f2Name = [SampleInfo.Pair2BamFile{iSample}];
  [p2Fwd, p2Rev, inBounds(3,iSample), firstPos, lastPos, inBounds(4,iSample)]...
      = createDepthChart(runtimeResultsPath, dataPath, refSeq,f2Name, cutoff_phred, firstPos, lastPos, qualOffset);

  rawDepth(iSample).P1Fwd = p1Fwd;
  rawDepth(iSample).P1Rev = p1Rev;
  rawDepth(iSample).P2Fwd = p2Fwd;
  rawDepth(iSample).P2Rev = p2Rev;
  clear p1Fwd p1Rev p2Fwd p2Rev
end
fprintf('\nDone.\n');

% Define region of interest to be analyzed
start = textscan(ROI, '%d:%d');
pos1= max(firstPos,start{1});
pos2 = min(lastPos,start{2});
goodPosIdx = pos1:pos2;
goodPosInd = false(1,length(refSeq)); goodPosInd(goodPosIdx) =true;
nPos = length(goodPosIdx);

% Align the sample sequence to the good bases of the alignment reference
refSeq(~goodPosInd) = [];

% Filter reads for good positions
for iSample = 1:nSample
  good = goodPosInd(inBounds(1,iSample):inBounds(2,iSample));
  rawDepth(iSample).P1Fwd(~good,:) = [];
  rawDepth(iSample).P1Rev(~good,:) = [];
  
  good = goodPosInd(inBounds(3,iSample):inBounds(4,iSample));
  rawDepth(iSample).P2Fwd(~good,:) = [];
  rawDepth(iSample).P2Rev(~good,:) = [];
end

%Reformat reference sequence bases
refBase=zeros(length(refSeq),1);
for iBase=1:length(refSeq)
    refBase(iBase)=refSeq(iBase);
end 

ref = nominal(refBase,{'1','2','3','4'});
refBase = nominal(ref,{'A','C','G','T'});

% Compile reference and error depth summaries
[ refDepth, errDepth, sumDepth] = compileRawDepth(rawDepth, ref);

% Export read count data in mat file
save(fullfile(runtimeResultsPath,'rawSeqData.mat'), 'refDepth', 'errDepth', 'sumDepth',...
  'rawDepth');

% Export position annotation data in mat file
SeqAnnot = dataset({goodPosIdx','AlignRefPos'}, {ref, 'AlignRefBase'}, {refBase, 'AlignRefChar'});
save(fullfile(runtimeResultsPath,'seqAnnot.mat'),'SeqAnnot');
end