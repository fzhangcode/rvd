function [ r, n, varargout ] = combineRawReads( errDepth, refDepth, rawDepth, isRef, sDiff, alpha)
% COMBINERAWREADS Combines forward/reverse and pairs while filtering
% differential read counts.

% Combine two forward reads and reverse reads due to high correlation.
errFwd = (errDepth.P1Fwd+errDepth.P2Fwd);
errRev = (errDepth.P1Rev+errDepth.P2Rev);
refFwd = (refDepth.P1Fwd+refDepth.P2Fwd);
refRev = (refDepth.P1Rev+refDepth.P2Rev);

% Filter strand errors by error differential
[ r, n, vFwdInd, vRevInd ] = ...
  filter_strand_by_diff( alpha, sDiff, isRef, errFwd, refFwd, errRev, refRev );

[nPos, nSample] = size(r);

% Update errDepth
errDepth.P1Fwd(vFwdInd,:) = 0;
errDepth.P2Fwd(vFwdInd,:) = 0;
errDepth.P1Rev(vRevInd,:) = 0;
errDepth.P2Rev(vRevInd,:) = 0;

% Update refDepth
refDepth.P1Fwd(vFwdInd,:) = 0;
refDepth.P2Fwd(vFwdInd,:) = 0;
refDepth.P1Rev(vRevInd,:) = 0;
refDepth.P2Rev(vRevInd,:) = 0;

% Update rawDepth
filDepth = struct([]);
for iSample = 1:nSample
    P1Fwd = zeros(nPos, 4);
    P2Fwd = zeros(nPos,4);
    P1Rev = zeros(nPos,4);
    P2Rev = zeros(nPos,4);
    for iPos = 1:nPos
        for baseIdx = 1:4
            if (iscell(rawDepth(iSample).P1Fwd))
                P1Fwd(iPos, baseIdx) = str2num(rawDepth(iSample).P1Fwd{iPos, baseIdx});
            else
                P1Fwd(iPos, baseIdx) = rawDepth(iSample).P1Fwd(iPos, baseIdx);
            end
                
            if (iscell(rawDepth(iSample).P2Fwd))
                P2Fwd(iPos, baseIdx) = str2num(rawDepth(iSample).P2Fwd{iPos, baseIdx});
            else
                P2Fwd(iPos, baseIdx) = rawDepth(iSample).P2Fwd(iPos, baseIdx);
            end
                
            if (iscell(rawDepth(iSample).P1Rev))
                P1Rev(iPos, baseIdx) = str2num(rawDepth(iSample).P1Rev{iPos, baseIdx});
            else
                P1Rev(iPos, baseIdx) = rawDepth(iSample).P1Rev(iPos, baseIdx);
            end
            
            if (iscell(rawDepth(iSample).P2Rev))
                P2Rev(iPos, baseIdx) = str2num(rawDepth(iSample).P2Rev{iPos, baseIdx});
            else
                P2Rev(iPos, baseIdx) = rawDepth(iSample).P2Rev(iPos, baseIdx);
            end
        end
     end
     filDepth(iSample).P1Fwd = P1Fwd;
     filDepth(iSample).P1Rev = P1Rev;
     filDepth(iSample).P2Fwd = P2Fwd;
     filDepth(iSample).P2Rev = P2Rev;
     clear P1Fwd P1Rev P2Fwd P2Rev
end

for iSample = 1:nSample
    filDepth(iSample).P1Fwd(vFwdInd,:) = 0;
    filDepth(iSample).P2Fwd(vFwdInd,:) = 0;
    filDepth(iSample).P1Rev(vRevInd,:) = 0;
    filDepth(iSample).P2Rev(vRevInd,:) = 0;
end

varargout{1} = errDepth;
varargout{2} = refDepth;
varargout{3} = filDepth;
end

