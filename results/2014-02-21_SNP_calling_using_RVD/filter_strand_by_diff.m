function [ r, n, varargout] = ...
	filter_strand_by_diff( alpha, sDiff, isRef, errFwd, refFwd, errRev, refRev )
%FILTER_STRAND_BY_DIFF Filter strand reads based on fwd/rev differential.


% Function checking
nargoutchk(2, 4);

[nPos, nSamp] = size(errFwd);

% Find differential error rates on the forward and reverse
% strands using reference data.
[fwdDiffInd, revDiffInd] = ...
	find_strand_diff(alpha, sDiff, errFwd(:,isRef), refFwd(:,isRef), ...
	errRev(:,isRef), refRev(:,isRef));

% Set the offending strand depths to zero
errFwd(fwdDiffInd,:) = zeros(sum(fwdDiffInd),nSamp);
refFwd(fwdDiffInd,:) = zeros(sum(fwdDiffInd),nSamp);
errRev(revDiffInd,:) = zeros(sum(revDiffInd),nSamp);
refRev(revDiffInd,:) = zeros(sum(revDiffInd),nSamp);

% Combine the forward and reverse strand depths after censoring
% error-positions
r = errFwd + errRev;
n = r + refFwd + refRev;

% Set output arguments
varargout(1) = {fwdDiffInd};
varargout(2) = {revDiffInd};

end

function [fwdDiffInd, revDiffInd] = ...
	find_strand_diff(alpha, sDiff, errFwd, refFwd, errRev, refRev)
% FIND_STRAND_DIFF Find differential error rates on the forward and reverse
% strands using reference data.

[nPos, nSamp] = size(errFwd);

% Calculate the forward and reverse strand MLE error rate across reference
% plexes
errFwd =sum(errFwd,2); refFwd =sum(refFwd,2);
errRev =sum(errRev,2); refRev =sum(refRev,2);

thetaFwd = errFwd./(errFwd+refFwd);
thetaRev = errRev./(errRev+refRev);

% Calculate the chi2 statistic for the contingency table Strand x Err/Ref

% Rows = Reference / Error; Columns = Forward / Reverse
Otbl = zeros(2,2,nPos);
Otbl(1,1,:) = refFwd;
Otbl(1,2,:) = refRev;
Otbl(2,1,:) = errFwd;
Otbl(2,2,:) = errRev;

Otot = squeeze(sum(sum(Otbl)));

% Calculate the expected values for the contingency table
errMarg = squeeze(sum(Otbl,2));
dirMarg = squeeze(sum(Otbl,1));

Etbl = zeros(2,2,nPos);
for i = 1:nPos
	Etbl(:,:,i) = errMarg(:,i)*dirMarg(:,i)'./Otot(i);
end

x = sum(sum(((Otbl - Etbl).^2)./Etbl,1),2);
x = squeeze(x);

% Identify positions with a significant error rate differential >0.005
x_thresh = chi2inv(1-alpha,1);
fwdDiffInd = x > x_thresh & thetaFwd-thetaRev > sDiff;
revDiffInd = x > x_thresh & thetaRev-thetaFwd > sDiff;
end

