function T = computeGrowthRatePairs(outputDir, cutoffFrame, posIdxRange, minLength, minParentLength, minDaughterLength, nCores, ycutoff)

% computeGrowthRatePairs
% Build mother-daughter growth-rate pairs from trackedCells.mat across positions
%
% Input:
% - outputDir: string/char, where the expInfoObj is loaded from
% - cutoffFrame: exclude cells with lastFrame >= cutoffFrame
% - posIdxRange: numeric indices into expInfoObj.getPositionList()
% - minLength/minParentLength/minDaughterLength: 
% - nCores: last argument passed to getGrowthRates, number of parallell
% workers
% - ycutoff: keep cells with ycoord > ycutoff. Use [] to not use.
%
% Output
% - T: table with columns: Mother, Daughter, Position

    if nargin < 2; cutoffFrame = []; end
    if nargin < 3; posIdxRange = []; end
    if nargin < 4 || isempty(minLength); minLength = 10; end
    if nargin < 5 || isempty(minParentLength); minParentLength = 0; end
    if nargin < 6 || isempty(minDaughterLength); minDaughterLength = 0; end
    if nargin < 7 || isempty(nCores); nCores = 8; end
    if nargin < 8; ycutoff = []; end

    expInfoObj = loadExpInfo(outputDir);
    posListAll = expInfoObj.getPositionList();
    if ~isempty(posIdxRange)
        posList = posListAll(posIdxRange);
    else
        posList = posListAll;
    end

    motherAll = [];
    daughterAll = [];
    posAll = strings(0,1);

    for pIdx = 1:numel(posList)
        pos = posList{pIdx};

        trackedFile = fullfile(outputDir, pos, "trackedCells.mat");
        if ~isfile(trackedFile)
            continue;
        end

        mCells = Cell.MCell.loadMCells(trackedFile);
        if isempty(mCells)
            continue;
        end

        mCells = mCells(:);
        n = numel(mCells);

        lastFrames = [mCells.lastFrame]';

        % y coordinate at end of track
        ycoord = nan(n, 1);
        for i = 1:n
            c = mCells(i);
            ycoord(i) = c.boundingBoxes(2, end);
        end

        sel = true(n, 1);

        if ~isempty(cutoffFrame)
            sel = sel & (lastFrames < cutoffFrame);
        end

        if ~isempty(ycutoff)
            sel = sel & (ycoord > ycutoff);
        end

        selCells = mCells(sel);
        if isempty(selCells)
            continue;
        end

        % Growth rates for selected cells
        [grs, genIndex] = getGrowthRates(selCells, minLength, minParentLength, minDaughterLength, nCores);

        if isempty(grs)
            continue;
        end

        % Lookup by cell id
        grs2 = nan(n, 1);
        grs2([selCells(genIndex).id]) = grs;

        % Mother-daughter pairs
        pairs = nan(n, 2);
        for i = 1:numel(grs)
            c = selCells(genIndex(i));
            pc = c.parent;
            if ~isempty(pc)
                pairs(c.id, :) = [grs2(pc.id), grs(i)];
            end
        end

        good = all(isfinite(pairs), 2);
        pairs = pairs(good, :);

        if isempty(pairs)
            continue;
        end

        motherAll = [motherAll; pairs(:,1)];
        daughterAll = [daughterAll; pairs(:,2)];
        posAll = [posAll; repmat(string(pos), size(pairs,1), 1)];
    end

    T = table(motherAll, daughterAll, posAll, 'VariableNames', {'Mother','Daughter','Position'});
end
