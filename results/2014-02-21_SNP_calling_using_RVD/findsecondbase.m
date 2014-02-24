function [ secondBase , secondBasePrc] = findsecondbase( ReadDepth, alignRefBase )
% Identify the second most frequent base by read depth

nPos = length(alignRefBase);

alignRefIdx = uint8(alignRefBase)';
alignRefChar = char(alignRefBase)';
alignRefInd = false(nPos,4);
for i = 1:nPos, alignRefInd(i,alignRefIdx(i)) = true; end

ntDepth = ReadDepth.P1Fwd + ReadDepth.P1Rev + ...
  ReadDepth.P2Fwd + ReadDepth.P2Rev;

nonrefDepth = ntDepth; nonrefDepth(alignRefInd) = NaN;

[secondBaseCount, secondBaseIdx] = max(nonrefDepth,[],2);
secondBasePrc = secondBaseCount./sum(ntDepth,2);
secondBase = nominal(secondBaseIdx,{'A','C','G','T'});

end

