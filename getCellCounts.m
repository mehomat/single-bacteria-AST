function T = getCellCounts(switchFrame, windowSize, labels, timeOffsets, tablefilename, varargin)

% getCellCounts
% Cell count analysis using a defined window (trailing by default)
%
% Good cell definition: isBadCell == 0 OR isBadCell == 4
% Per frame requirement: badSegmentations(frameIdx) == 0
%
% Input arguments:
% - switchFrame: scalar, frame index of media switch
% - windowSize: scalar, number of frames in the window
% - labels: string array (2 elements)
% - timeOffsets: numeric row vector, frames after media switch
% - tablefilename: string/char, output CSV filename to write
%
% Name-value pairs:
% - 'UseTrailingWindow': logical, default true (false --> centered window)
% - 'CodeDir': string/char, default 'code' (kept for compatibility)
% - 'OutputDir': string/char, required

%% Parse inputs

p = inputParser;
p.addParameter('CodeDir', 'code');
p.addParameter('OutputDir', '');
p.addParameter('UseTrailingWindow', true, @(x)islogical(x) || isnumeric(x));
p.parse(varargin{:});

outputDir = p.Results.OutputDir;
useTrailingWindow = logical(p.Results.UseTrailingWindow);

% Ensure row vector
timeOffsets = timeOffsets(:).';

%% Load experiment data

expInfoObj = loadExpInfo(outputDir);
posList = expInfoObj.getPositionList();
param = expInfoObj.getParameters();
nGrowthChannels = param.nGrowthChannels - length(param.emptyChannel);

%% Experimental settings

strainNames = [labels(1), labels(2)];
timePoints = switchFrame + timeOffsets;

% Define frame range that actually needs to be computed
if useTrailingWindow
    minNeeded = min(timePoints) - (windowSize - 1);
    maxNeeded = max(timePoints);
else
    halfw = floor(windowSize/2);
    minNeeded = min(timePoints) - halfw;
    maxNeeded = max(timePoints) + halfw;
end
minNeeded = max(1, minNeeded);

%% Initialize storage

nPos = numel(posList);
allCellCounts  = cell(nPos, 1);
allStrains = cell(nPos, 1);
allTrapNumbers = cell(nPos, 1);

%% Loop over positions

for pIdx = 1:nPos

    posName = posList{pIdx};
    trackedPath = fullfile(outputDir, posName, 'trackedCells.mat');

    % Load tracked cells
    try
        [mCells, ~] = Cell.MCell.loadMCells(trackedPath);
    catch
        fprintf('Could not load mCells for %s\n', posName);
        allCellCounts{pIdx}  = [];
        allTrapNumbers{pIdx} = [];
        allStrains{pIdx} = [];
        continue;
    end

    if isempty(mCells)
        allCellCounts{pIdx} = [];
        allTrapNumbers{pIdx} = [];
        allStrains{pIdx} = [];
        continue;
    end

    % Precompute per cell arrays
    gc = [mCells.growthChannel];
    bf = [mCells.birthFrame];
    lf = [mCells.lastFrame];
    bad = [mCells.isBadCell];

    % Used for good/bad cell
    goodGlobal = (bad == 0) | (bad == 4);

    % Determine max frame one can access
    maxFromLife = max(lf);
    maxF = min(maxNeeded, maxFromLife);

    % Ensure switchFrame and requested frames are within computed range
    if switchFrame > maxF

        % No usable data for requested timepoints
        allCellCounts{pIdx}  = [];
        allTrapNumbers{pIdx} = [];
        allStrains{pIdx}     = [];
        continue;
    end

    % counts(trap, frame) for frames 1..maxF
    counts = zeros(nGrowthChannels, maxF, 'single');

    % Build counts matrix once
    for ci = 1:numel(mCells)
        if ~goodGlobal(ci)
            continue;
        end

        t = gc(ci);
        if t < 1 || t > nGrowthChannels
            continue;
        end

        % Segmentation vector coverage: frames are bf(ci) ... bf(ci)+numel(seg)-1
        seg = mCells(ci).badSegmentations;
        if isempty(seg)
            continue;
        end
        segMaxFrame = bf(ci) + numel(seg) - 1;

        % Frame interval where this cell could contribute
        f1 = max([bf(ci), minNeeded, 1]);
        f2 = min([lf(ci), maxF, segMaxFrame]);

        if f2 < f1
            continue;
        end

        idx1 = f1 - bf(ci) + 1;
        idx2 = f2 - bf(ci) + 1;

        ok = (seg(idx1:idx2) == 0);
        if ~any(ok)
            continue;
        end

        frames = f1:f2;
        frames = frames(ok);

        % Increment counts for those frames
        counts(t, frames) = counts(t, frames) + 1;
    end

    % Identify valid traps at switchFrame (>= 3 good, well segmented cells at switch)
    trapCountsAtSwitch = counts(:, switchFrame).';
    validTraps = find(trapCountsAtSwitch >= 3);
    nValidTraps = numel(validTraps);

    if nValidTraps == 0
        allCellCounts{pIdx}  = [];
        allTrapNumbers{pIdx} = [];
        allStrains{pIdx} = [];
        continue;
    end

    % Prepare output arrays
    trapCounts  = nan(nValidTraps, numel(timePoints));
    trapNumbers = (pIdx - 1) * nGrowthChannels + validTraps(:);

    % Compute windowed means using cumulative sums
    for ii = 1:nValidTraps
        tID = validTraps(ii);

        x = double(counts(tID, :)); % 1 x maxF
        cs = cumsum([0, x]); % 1 x (maxF+1)

        for jj = 1:numel(timePoints)
            f = timePoints(jj);

            if f < 1 || f > maxF
                trapCounts(ii, jj) = NaN;
                continue;
            end

            if useTrailingWindow
                fStart = max(f - windowSize + 1, 1);
                fEnd = f;
            else
                halfw = floor(windowSize/2);
                fStart = max(f - halfw, 1);
                fEnd = min(f + halfw, maxF);
            end

            % Mean in [fStart, fEnd]
            sumWin = cs(fEnd + 1) - cs(fStart);
            nWin = (fEnd - fStart + 1);
            trapCounts(ii, jj) = sumWin / nWin;
        end
    end

    % Assign strain label by position (first half vs second half)
    if pIdx <= nPos/2
        strain = strainNames(1);
    else
        strain = strainNames(2);
    end

    allCellCounts{pIdx}  = trapCounts;
    allTrapNumbers{pIdx} = trapNumbers;
    allStrains{pIdx}     = repmat(strain, size(trapCounts, 1), 1);
end

%% Write table
if all(cellfun(@isempty, allCellCounts))
    warning('No traps aggregated. Returning empty table.');
    T = table();
    if strlength(string(tablefilename)) > 0
        writetable(T, tablefilename);
    end
    return;
end

allCellCountsMat = cell2mat(allCellCounts);
allTrapNumbersCol = vertcat(allTrapNumbers{:});
allStrainsCol = vertcat(allStrains{:});

varNames = strcat(string(timeOffsets), " min");
T = array2table(allCellCountsMat, 'VariableNames', varNames);
T.TrapNumber = allTrapNumbersCol;
T.Strain = allStrainsCol;

writetable(T, tablefilename);
fprintf('Wrote %s with %d traps x %d timepoints.\n', tablefilename, size(T,1), numel(timeOffsets));
end
