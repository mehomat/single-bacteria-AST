function T = getSingleLineageCellAreasFig1(pIdx, trapIdx, endFrame, varargin)

% T = getSingleLineageCellAreasFig1(pIdx, trapIdx, endFrame, 'OutputDir', outputDir)
% One lineage (lower child) from trap mother, using the same trap-mother logic
% as getSingleLineageCellAreas. Returns table(Frame, Area).

p = inputParser;
p.addParameter('OutputDir', '');
p.parse(varargin{:});
outputDir = p.Results.OutputDir;

expInfoObj = loadExpInfo(outputDir);
posList = expInfoObj.getPositionList();

trackedPath = fullfile(outputDir, posList{pIdx}, 'trackedCells.mat');
mCells = Cell.MCell.loadMCells(trackedPath);

trapMotherCellIds = find_trap_mother_ids_local(mCells, trapIdx);
cid = trapMotherCellIds(1);
c = mCells(cid);

frames = [];
areas  = [];

while ~isempty(c) && (c.isBadCell == 0 || c.isBadCell == 2) && c.birthFrame <= endFrame

    f1 = c.birthFrame;
    f2 = min(c.lastFrame, endFrame);
    n  = f2 - f1 + 1;

    cframes = f1:f2;
    careas = c.areas(1:n).';
    bad = c.badSegmentations(1:n).';
    careas(bad > 0) = NaN;

    frames = [frames, cframes];
    areas = [areas, careas.'];

    if isempty(c.descendants)
        break;
    end

    if numel(c.descendants) < 2
        c = c.descendants(1);
    else
        dc1 = c.descendants(1);
        dc2 = c.descendants(2);

        idx1 = find(dc1.badSegmentations == 0, 1);
        idx2 = find(dc2.badSegmentations == 0, 1);

        dc1_y = NaN; if ~isempty(idx1), dc1_y = dc1.centroids(idx1,2); end
        dc2_y = NaN; if ~isempty(idx2), dc2_y = dc2.centroids(idx2,2); end

        if dc1_y > dc2_y
            c = dc1; % lower child
        else
            c = dc2;
        end
    end
end

T = table(frames(:), areas(:), 'VariableNames', ["Frame","Area"]);

end

%% Local helper
function trapMotherCellIds = find_trap_mother_ids_local(mCells, trapIdx)

trapMotherCellIds = [];

birthFrames = [mCells.birthFrame];
lastFrames  = [mCells.lastFrame];
cellTraps = [mCells.growthChannel];

indTrap = (cellTraps == trapIdx);

cellYcoord = nan(size(cellTraps));
cellLengths = nan(size(cellTraps));
shifty = nan(size(cellTraps));

for i = 1:numel(mCells)
    f = find(mCells(i).badSegmentations == 0, 1);
    if ~isempty(f)
        bb = mCells(i).boundingBoxes(:, f);
        cellYcoord(i)  = bb(2);
        cellLengths(i) = bb(4);

        goodFrames = mCells(i).badSegmentations == 0;
        if nnz(goodFrames) > 1
            shifty(i) = mean(diff(mCells(i).centroids(goodFrames, 2)), 'omitnan');
        end
    end
end

% dy logic
lengthsPool = cellLengths(cellLengths > 0);
if isempty(lengthsPool) || all(isnan(lengthsPool))
    dy = 0;
else
    dy = quantile(lengthsPool, 0.1) / 2;
end

meanshifty = mean(shifty, 'omitnan');
if isnan(meanshifty)
    meanshifty = 0;
end

if meanshifty < 0
    cellYcoord = cellYcoord + cellLengths;
end

if meanshifty > 0
    yMin = min(cellYcoord(indTrap), [], 'omitnan');
    yCutOff = yMin + dy;
    trapMotherCellIds = find(indTrap & cellYcoord < yCutOff);
else
    yMax = max(cellYcoord(indTrap), [], 'omitnan');
    yCutOff = yMax - dy;
    trapMotherCellIds = find(indTrap & cellYcoord > yCutOff);
end

[~, sortInd] = sort(birthFrames(trapMotherCellIds));
trapMotherCellIds = trapMotherCellIds(sortInd);

end


