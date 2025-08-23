function clusters = findClusters(B, conn)
    if nargin < 2, conn = 4; end
    CC = bwconncomp(B ~= 0, conn);
    clusters = CC.PixelIdxList;
end
