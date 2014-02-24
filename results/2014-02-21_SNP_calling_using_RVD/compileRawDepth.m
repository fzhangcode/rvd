function [ refDepth, errDepth, sumDepth ] = compileRawDepth(rawDepth, ref)
%COMPILERAWDEPTH Converts raw depths to reference and error depth counts.


nSample = length(rawDepth);
nPos = size(rawDepth(1).P1Fwd,1);

refBaseIdx = uint8(ref);
refDepth.P1Fwd = zeros(nPos,nSample); refDepth.P1Rev = zeros(nPos,nSample);
refDepth.P2Fwd = zeros(nPos,nSample); refDepth.P2Rev = zeros(nPos,nSample);

errDepth.P1Fwd = zeros(nPos,nSample); errDepth.P1Rev = zeros(nPos,nSample);
errDepth.P2Fwd = zeros(nPos,nSample); errDepth.P2Rev = zeros(nPos,nSample);

sumDepth.P1Fwd = zeros(nPos,nSample); sumDepth.P1Rev = zeros(nPos,nSample);
sumDepth.P2Fwd = zeros(nPos,nSample); sumDepth.P2Rev = zeros(nPos,nSample);

for iSample = 1:nSample
  for iPos = 1:nPos
    refInd = false(1,4); refInd(refBaseIdx(iPos)) = true;
    
    if(iscell(rawDepth(iSample).P1Fwd(iPos,refInd)))
        refDepth.P1Fwd(iPos, iSample) = str2num(rawDepth(iSample).P1Fwd{iPos, refInd});
        errDepth.P1Fwd(iPos, iSample) = sum(cellfun(@str2num,rawDepth(iSample).P1Fwd(iPos, ~refInd)));
    else
        refDepth.P1Fwd(iPos, iSample) = rawDepth(iSample).P1Fwd(iPos, refInd);
        errDepth.P1Fwd(iPos, iSample) = sum(rawDepth(iSample).P1Fwd(iPos, ~refInd));
    end
    
    if(iscell(rawDepth(iSample).P2Fwd(iPos,refInd)))
        refDepth.P2Fwd(iPos, iSample) = str2num(rawDepth(iSample).P2Fwd{iPos, refInd});
        errDepth.P2Fwd(iPos, iSample) = sum(cellfun(@str2num,rawDepth(iSample).P2Fwd(iPos, ~refInd)));
    else
        refDepth.P2Fwd(iPos, iSample) = rawDepth(iSample).P2Fwd(iPos, refInd);
        errDepth.P2Fwd(iPos, iSample) = sum(rawDepth(iSample).P2Fwd(iPos, ~refInd));
    end
    
    if(iscell(rawDepth(iSample).P1Rev(iPos,refInd)))
        refDepth.P1Rev(iPos, iSample) = str2num(rawDepth(iSample).P1Rev{iPos, refInd});
        errDepth.P1Rev(iPos, iSample) = sum(cellfun(@str2num,rawDepth(iSample).P1Rev(iPos, ~refInd)));
    else
        refDepth.P1Rev(iPos, iSample) = rawDepth(iSample).P1Rev(iPos, refInd);
        errDepth.P1Rev(iPos, iSample) = sum(rawDepth(iSample).P1Rev(iPos, ~refInd));
    end
    
    if(iscell(rawDepth(iSample).P2Rev(iPos,refInd)))
        refDepth.P2Rev(iPos, iSample) = str2num(rawDepth(iSample).P2Rev{iPos, refInd});
        errDepth.P2Rev(iPos, iSample) = sum(cellfun(@str2num,rawDepth(iSample).P2Rev(iPos, ~refInd)));
    else
        refDepth.P2Rev(iPos, iSample) = rawDepth(iSample).P2Rev(iPos, refInd);
        errDepth.P2Rev(iPos, iSample) = sum(rawDepth(iSample).P2Rev(iPos, ~refInd));
    end
    
  end
  sumDepth.P1Fwd(:, iSample) = refDepth.P1Fwd(:,iSample) + errDepth.P1Fwd(:,iSample);
  sumDepth.P1Rev(:, iSample) = refDepth.P1Rev(:,iSample) + errDepth.P1Rev(:,iSample);
  sumDepth.P2Fwd(:, iSample) = refDepth.P2Fwd(:,iSample) + errDepth.P2Fwd(:,iSample);
  sumDepth.P2Rev(:, iSample) = refDepth.P2Rev(:,iSample) + errDepth.P2Rev(:,iSample);
end


end

