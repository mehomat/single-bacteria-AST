function T = getSingleLineageCellAreas(labels, positionRes, trapRes, positionSus, trapSus, lineagefilename, outputDir, varargin)

% getSingleLineageCellAreas
% Stores the cell areas for two cell lineages, one resistant and one susceptible
% Lineage choice is always the lower child in the trap
%
% Label handling:
% - labels must contain exactly 2 entries 
% - One label must end with 'resistant' and the other must end with 'susceptible'
% - The function assigns labelRes/labelSus from that suffix
% - Uses (positionRes, trapRes) for the resistant lineage and
%  (positionSus,trapSus) for the susceptible lineage 
%
% Inputs:
% - labels: string array with 2 entries
% - positionRes, trapRes: scalars
% - positionSus, trapSus: scalars
% - lineagefilename: output CSV filename
% - outputDir: string/char, where the expInfoObj is loaded from
%
% Name-value pairs:
% - 'SwitchFrame': scalar, default Inf
% - 'YThresh': scalar or [], default []
%
% Output:
% - T: table with the cell areas for the two cell lineages

%% Parse inputs
p = inputParser;
p.addParameter('SwitchFrame', Inf);
p.addParameter('YThresh', []);
p.parse(varargin{:});

switchFrame = p.Results.SwitchFrame;
yThresh = p.Results.YThresh;

%% Parse labels

labels = strtrim(labels); %T rim
[labelRes, labelSus] = parseResSusLabels(labels);

%% Load experiment data

expInfoObj = loadExpInfo(outputDir);
posList = expInfoObj.getPositionList();

%% Helper: Growth rate estimate before switch

    function gr = getGrowthRateBeforeSwitch_local(cellObj, swFrame)

        gr = NaN;

        if isempty(cellObj) || isempty(cellObj.boundingBoxes) || isempty(cellObj.badSegmentations)
            return;
        end

        nFramesCell = numel(cellObj.badSegmentations);
        if isfinite(swFrame)
            maxF = min(nFramesCell, max(1, floor(swFrame - cellObj.birthFrame + 1)));
        else
            maxF = nFramesCell;
        end

        good = (cellObj.badSegmentations(1:maxF) == 0);
        if nnz(good) < 3
            return;
        end

        L = double(cellObj.boundingBoxes(4, 1:maxF)).';
        fIdx = double((1:maxF).');

        keep = good(:) & isfinite(L) & (L > 0);
        if nnz(keep) < 3
            return;
        end

        pp = polyfit(fIdx(keep), log(L(keep)), 1);
        gr = pp(1);
    end

%% Helper: find trap mother IDs

    function trapMotherCellIds = find_trap_mother_ids_local(mCells, trapIdx, swFrame, yThr)

        if nargin < 4
            yThr = [];
        end

        trapMotherCellIds = [];

        if isempty(mCells)
            return;
        end

        birthFrames = [mCells.birthFrame];
        lastFrames = [mCells.lastFrame];
        cellTraps = [mCells.growthChannel];

        indTrap = (cellTraps == trapIdx);
        if ~any(indTrap)
            return;
        end

        cellYcoord  = nan(size(cellTraps));
        cellLengths = nan(size(cellTraps));
        shifty = nan(size(cellTraps));
        cellGR = nan(size(cellTraps));

        for i = 1:numel(mCells)
            f = find(mCells(i).badSegmentations == 0, 1); % first good frame
            if ~isempty(f)
                bb = mCells(i).boundingBoxes(:, f);
                cellYcoord(i)  = bb(2);
                cellLengths(i) = bb(4);

                goodFrames = mCells(i).badSegmentations == 0;
                if nnz(goodFrames) > 1
                    shifty(i) = mean(diff(mCells(i).centroids(goodFrames, 2)), 'omitnan');
                end
            end

            cellGR(i) = getGrowthRateBeforeSwitch_local(mCells(i), swFrame);

        end

        lengthsPool = cellLengths(birthFrames < swFrame & cellLengths > 0);
        if isempty(lengthsPool) || all(isnan(lengthsPool))
            dy = 0;
        else
            dy = quantile(lengthsPool, 0.1) / 2;
        end

        meanshifty = mean(shifty(lastFrames < swFrame), 'omitnan');
        if isnan(meanshifty)
            meanshifty = 0;
        end

        if meanshifty < 0
            cellYcoord = cellYcoord + cellLengths;
        end

        if ~isempty(yThr)
            if meanshifty < 0
                cellYcoord(cellYcoord > yThr) = NaN;
            else
                cellYcoord(cellYcoord < yThr) = NaN;
            end
        end

        if meanshifty > 0
            yMin = min(cellYcoord(indTrap), [], 'omitnan');
            if isnan(yMin), return; end
            yCutOff = yMin + dy;
            trapMotherCellIds = find(indTrap & cellYcoord < yCutOff);
        else
            yMax = max(cellYcoord(indTrap), [], 'omitnan');
            if isnan(yMax), return; end
            yCutOff = yMax - dy;
            trapMotherCellIds = find(indTrap & cellYcoord > yCutOff);
        end

        if isempty(trapMotherCellIds)
            return;
        end

        [~, sortInd] = sort(birthFrames(trapMotherCellIds));
        trapMotherCellIds = trapMotherCellIds(sortInd);
    end

%% Helper: extract a single lineage (lower child) starting from trap mother
    
    function [frames, areas, motherId] = extract_lineage(pIdx, trapIdx)

        trackedPath = fullfile(outputDir, posList{pIdx}, 'trackedCells.mat');
        mCells = Cell.MCell.loadMCells(trackedPath);

        frames = [];
        areas = [];
        motherId = NaN;

        if isempty(mCells)
            return;
        end

        trapMotherCellIds = find_trap_mother_ids_local(mCells, trapIdx, switchFrame, yThresh);
        if isempty(trapMotherCellIds)
            return;
        end

        cid = trapMotherCellIds(1);
        c = mCells(cid);
        motherId = c.id;

        while ~isempty(c) && (c.isBadCell == 0 || c.isBadCell == 2)
            cframes = c.birthFrame:c.lastFrame;

            careas = c.areas(:).';
            bad = c.badSegmentations(:).';
            careas(bad > 0) = NaN;

            frames = [frames, cframes];
            areas = [areas,  careas];

            if isempty(c.descendants)
                c = [];
            else
                if numel(c.descendants) < 2
                    c = c.descendants(1);
                else
                    dc1 = c.descendants(1);
                    idx1 = find(dc1.badSegmentations == 0, 1);
                    dc1_y = NaN; if ~isempty(idx1), dc1_y = dc1.centroids(idx1,2); end

                    dc2 = c.descendants(2);
                    idx2 = find(dc2.badSegmentations == 0, 1);
                    dc2_y = NaN; if ~isempty(idx2), dc2_y = dc2.centroids(idx2,2); end

                    if all(isnan([dc1_y, dc2_y]))
                        break;
                    elseif dc1_y > dc2_y
                        c = dc1;
                    else
                        c = dc2;
                    end
                end
            end
        end
    end

%% Extract lineages 

[framesRes, areasRes, motherRes] = extract_lineage(positionRes, trapRes);
[framesSus, areasSus, motherSus] = extract_lineage(positionSus, trapSus);

%% Combine and save

Tres = table(framesRes(:), areasRes(:), repmat(labelRes, numel(framesRes), 1), ...
    'VariableNames', ["Frame","Area","Strain"]);
Tsus = table(framesSus(:), areasSus(:), repmat(labelSus, numel(framesSus), 1), ...
    'VariableNames', ["Frame","Area","Strain"]);

% Add columns
Tres.Position = repmat(positionRes, height(Tres), 1);
Tres.Trap = repmat(trapRes, height(Tres), 1);

Tsus.Position = repmat(positionSus, height(Tsus), 1);
Tsus.Trap = repmat(trapSus, height(Tsus), 1);

T = [Tres; Tsus];
writetable(T, lineagefilename);

end

%% Helper

function [labelRes, labelSus] = parseResSusLabels(labels)

    low1 = lower(labels(1));
    low2 = lower(labels(2));

    isRes1 = endsWith(low1, "resistant");
    isSus1 = endsWith(low1, "susceptible");
    isRes2 = endsWith(low2, "resistant");
    isSus2 = endsWith(low2, "susceptible");

    if isRes1
        labelRes = canonicalSuffix(labels(1), "resistant");
        labelSus = canonicalSuffix(labels(2), "susceptible");
    else
        labelRes = canonicalSuffix(labels(2), "resistant");
        labelSus = canonicalSuffix(labels(1), "susceptible");
    end
end

function lbl = canonicalSuffix(lbl, suffix)

    expr = "\s*" + suffix + "\s*$";
    lbl = regexprep(lbl, expr, " " + suffix, "ignorecase");

end




