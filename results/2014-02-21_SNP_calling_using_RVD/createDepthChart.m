
function [seqDepthFwd, seqDepthRev, sample_Start, firstPos, lastPos, sample_Stop]...
    = createDepthChart(runtimeResultsPath, dataPath, refSeq, fName, cutoff_phred, firstPos, lastPos,qualOffset)
% Creates Depth Chart from pileup
%
% Inputs:
%   runtimeResultsPath - directory to store temporary data
%   dataPath - path to data directory
%   refSeq - reference sequence
%   fName - name of bam sequence sample file
%   cutoff_phred - base quality threshold for filtering low quality reads
%   firstPos - first postion with > 1000 depth in all the previously reviewed BAM files
%   lastPos - last position with > 1000 depth in all the previously reviewed BAM files
%
% Outputs:
%   seqDepthFwd
%   seqDepthRev
%   sample_Start - first in bounds (>1000 depth) position in current BAM file
%   firstPos - first position that so far has >1000 depth in all reviewed
%       BAM files (including the current)
%   lastPos - last position that so far has >1000 depth in all reviewed
%       BAM files (including the current)
%   sample_Stop - last in bounds (>1000 depth) position in current BAM file
%
%
% Author: Anna Cushing
% Created: 
% Modified: November 7, 2012

%c\Converts BAM file to pileup and then to depth chart

sample_Start = firstPos;
sample_Stop = lastPos;

%needs to be sorted bam file
pileupName = strcat(fName, '.pileup');

filename = fullfile(dataPath, fName);

% Check that the file exists
if ~exist(filename, 'file'), error(['File ',filename,' does not exist.']); end

sortedName = fullfile(runtimeResultsPath, sprintf('%s.st', fName));
[status,result]=system(sprintf('samtools sort "%s" "%s"', filename, sortedName));
assert(status == 0, result);
sortedName = sprintf('"%s".bam',sortedName);
[status,result]=system(sprintf('samtools mpileup -d 10000000 "%s" > "%s"', sortedName,fullfile(runtimeResultsPath,pileupName)));

%Check that pileup conversion was successful
assert(status == 0, result);

%Make basic DC.
fid=fopen(fullfile(runtimeResultsPath, pileupName));
out_filename=[fullfile(runtimeResultsPath, fName) '.DC.txt'];

 %get length of pileup file and break it into blocks for reading
 [s,w] = system(sprintf('wc -l "%s"', fullfile(runtimeResultsPath,pileupName)));
 L=textscan(w,'%n',1);
 nblocks=ceil(L{1}/1000);

 %open out file to write depth chart to
 fid_out=fopen(out_filename, 'wt');
 
 %write header for depth chart
 fprintf(fid_out,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',...
'pos','refb','A_All','C_All','G_All','T_All','N_All','Tot_All','A_F','C_F','G_F','T_F','N_F','Tof_F','A_R','C_R','G_R','T_R','N_R','Tot_R');
 
 %read in blocks and scan out chr, pos, refb, depth, and reads
 out_count=zeros(L{1},10);
 for block=1:nblocks
    C=textscan(fid, '%s %s %s %n %s %s',1000, 'delimiter','\t','Bufsize',10000000);
    chr=C{1};pos=C{2};refb=C{3}; depth=C{4}; read=C{5}; qc=C{6}; %clear C
    if block==1
        sample_Start = pos{1};
        while depth(sample_Start) < 1000
            sample_Start = sample_Start +1;
        end
        sample_Stop = sample_Start;
    end
    
    %ignore insertion/deletion reads
    read=regexprep(read,'\^.|\$|-(\d)+(??[ACTGNacgtn]{$1})|\+(\d)+(??[ACTGNacgtn]{$1})','');
    c=cellfun(@length,read);
    l=length(find((depth-c)~=0));
    bases='ACGTNacgtn';
 
    block_start = (block-1)*1000+1;
    block_end = min(block*1000,L{1});
    block_len = min(block*1000,L{1}) - (block-1)*1000;
    
    for m=1:block_len
       if depth(m) > 1000
           sample_Stop = pos{m};
       end
       QC=char(qc{m});
       QC=double(QC); QC=QC-qualOffset;
       read{m}=read{m}(QC>cutoff_phred);
    end
    
    for n=1:10
     r=regexp(read,bases(n));
     out_count(block_start:block_end,n)=cellfun(@length,r);
    end
 end
 sample_Start = str2double(sample_Start);
 sample_Stop = str2double(sample_Stop);
 A_F = out_count(sample_Start:sample_Stop,1);
 C_F = out_count(sample_Start:sample_Stop,2);
 G_F = out_count(sample_Start:sample_Stop,3);
 T_F = out_count(sample_Start:sample_Stop,4);
 N_F = out_count(sample_Start:sample_Stop,5);
 
 A_R = out_count(sample_Start:sample_Stop,6);
 C_R = out_count(sample_Start:sample_Stop,7);
 G_R = out_count(sample_Start:sample_Stop,8);
 T_R = out_count(sample_Start:sample_Stop,9);
 N_R = out_count(sample_Start:sample_Stop,10);
 
 
 tot_F = A_F + C_F + G_F + T_F + N_F;
 tot_R  = A_R + C_R + G_R + T_R + N_R;
 tot_All  = tot_F + tot_R;
 A_All  = A_F + A_R;
 C_All  = C_F + C_R;
 G_All  = G_F + G_R;
 T_All  = T_F + T_R;
 N_All  = N_F + N_R;
 
 
 seqDepthFwd = [A_F C_F G_F T_F tot_F];
 seqDepthRev = [A_R C_R G_R T_R tot_R];

 %print ouput to file
  firstPos = max(sample_Start,firstPos);
  lastPos = min(sample_Stop,lastPos);
  
  for m=sample_Start:sample_Stop
      fprintf(fid_out, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t',num2str(m),refSeq(m), num2str(A_All(m)), num2str(C_All(m)),num2str(G_All(m)),num2str(T_All(m)),num2str(N_All(m)),num2str(tot_All(m)));
      for n=1:5
        fprintf(fid_out,'%s\t',num2str(out_count(m,n)));
      end
      fprintf(fid_out,'%s\t', num2str(tot_F(m)));
      for n=6:10
        fprintf(fid_out,'%s\t',num2str(out_count(m,n)));
      end
      fprintf(fid_out,'%s\t',num2str(tot_R(m)));
      fprintf(fid_out, '\n');
 end
 fclose(fid_out);
end
 