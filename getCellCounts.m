function T = getCellCounts(switchFrame, windowSize, labels, timeOffsets, tablefilename, posRange, outputDir, varargin)

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
% - posRange: 1x2 cell array where:
%  posRange{1} = numeric vector of positions for labels(1)
%  posRange{2} = numeric vector of positions for labels(2)
% - outputDir: string/char, where the expInfoObj is loaded from
%
% Name-value pairs:
% - 'UseTrailingWindow': logical, default true (false --> centered window)
%
% Output:
% - T: table of cell counts

%% Parse inputs

p = inputParser;
p.addParameter('UseTrailingWindow', true, @(x)islogical(x) || isnumeric(x));
p.parse(varargin{:});

useTrailingWindow = logical(p.Results.UseTrailingWindow);

% Ensure row vector
timeOffsets = timeOffsets(:).';

%% Load experiment data

expInfoObj = loadExpInfo(outputDir);
posList = expInfoObj.getPositionList();
param = expInfoObj.getParameters();
nGrowthChannels = param.nGrowthChannels - length(param.emptyChannel);

%% Experimental settings

timePoints = switchFrame + timeOffsets;

% Strain names
strainNames = [labels(1), labels(2)];

% Position ranges
posLabel1 = unique(posRange{1}(:).');
posLabel2 = unique(posRange{2}(:).');

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

% Store trap inclusion counts per position
nNonEmptyTraps_perPos = nan(nPos, 1);
nIncludedTraps_perPos = nan(nPos, 1);

%% Loop over positions

for pIdx = 1:nPos

    mCells = expInfoObj.getMCells(pIdx);

    if isempty(mCells)
        continue;
    end

    % Get valid traps
    S = getValidTraps(expInfoObj, pIdx, switchFrame, 1, ...
        'Traps', 1:nGrowthChannels);

    if ~isempty(S)
        nNonEmptyTraps_perPos(pIdx) = S(1).nNonEmptyTraps;
        nIncludedTraps_perPos(pIdx) = S(1).nIncludedTraps;
    else
        nNonEmptyTraps_perPos(pIdx) = 0;
        nIncludedTraps_perPos(pIdx) = 0;
    end
    
    validTraps = [S([S.isValid]).trap]';

    if isempty(validTraps)
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
        continue;
    end

    % counts(trap, frame) for frames 1..maxF
    counts = zeros(nGrowthChannels, maxF, 'single');

    % Build counts matrix
    for ci = 1:numel(mCells)
        if ~goodGlobal(ci)
            continue;
        end

        t = gc(ci);

        % Segmentation vector coverage: frames are bf(ci) ... bf(ci)+numel(seg)-1
        seg = mCells(ci).badSegmentations;

        % Frame interval where this cell could contribute
        f1 = max(bf(ci), minNeeded);
        f2 = min(lf(ci), maxF);

        if f2 < f1
            continue;
        end

        idx1 = f1 - bf(ci) + 1;
        idx2 = f2 - bf(ci) + 1;

        ok = (seg(idx1:idx2) == 0);
        frames = f1:f2;
        frames = frames(ok);

        % Increment counts for those frames
        counts(t, frames) = counts(t, frames) + 1;

    end

    nValidTraps = numel(validTraps);

    if nValidTraps == 0
        continue;
    end

    % Prepare output arrays
    trapCounts = nan(nValidTraps, numel(timePoints));
    trapNumbers = (pIdx - 1) * nGrowthChannels + validTraps(:);

    % Compute windowed means using cumulative sums
    for ii = 1:nValidTraps
        tID = validTraps(ii);

        x = double(counts(tID, :)); % 1 x maxF
        cs = cumsum([0, x]); % 1 x (maxF+1)

        for jj = 1:numel(timePoints)
            f = timePoints(jj);

            if f > maxF
                trapCounts(ii, jj) = NaN;
                continue;
            end

            if useTrailingWindow
                fStart = max(f - windowSize + 1, 1);
                fEnd = f;
            else
                fStart = max(f - halfw, 1);
                fEnd = min(f + halfw, maxF);
            end

            % Mean in [fStart, fEnd]
            sumWin = cs(fEnd + 1) - cs(fStart);
            nWin = (fEnd - fStart + 1);
            trapCounts(ii, jj) = sumWin / nWin;
        end
    end

    % Assign strain by position ranges
    if ismember(pIdx, posLabel1)
        strain = strainNames(1); % labels(1)
    elseif ismember(pIdx, posLabel2)
        strain = strainNames(2); % labels(2)
    else
        strain = ""; % or unassigned
    end

    allCellCounts{pIdx} = trapCounts;
    allTrapNumbers{pIdx} = trapNumbers;
    allStrains{pIdx} = repmat(strain, size(trapCounts, 1), 1);
end

%%

% Totals per strain over the requested posRange 
idxLabel1 = ismember((1:length(posList))', posLabel1);
idxLabel2 = ismember((1:length(posList))', posLabel2);

totalNonEmpty_1 = sum(nNonEmptyTraps_perPos(idxLabel1), 'omitnan');
totalIncluded_1 = sum(nIncludedTraps_perPos(idxLabel1), 'omitnan');

totalNonEmpty_2 = sum(nNonEmptyTraps_perPos(idxLabel2), 'omitnan');
totalIncluded_2 = sum(nIncludedTraps_perPos(idxLabel2), 'omitnan');

if totalNonEmpty_1 == 0
    pctIncluded_1 = NaN;
else
    pctIncluded_1 = 100 * totalIncluded_1 / totalNonEmpty_1;
end

if totalNonEmpty_2 == 0
    pctIncluded_2 = NaN;
else
    pctIncluded_2 = 100 * totalIncluded_2 / totalNonEmpty_2;
end

% Build + write summary CSV

Summary = table( ...
    strainNames(:), ...
    [totalNonEmpty_1; totalNonEmpty_2], ...
    [totalIncluded_1; totalIncluded_2], ...
    [pctIncluded_1; pctIncluded_2], ...
    'VariableNames', {'Strain','TotalNonEmptyTraps','TotalIncludedTraps','PercentIncludedTraps'});

summaryFilename = replace(string(tablefilename), ".csv", "_Summary.csv");
writetable(Summary, summaryFilename);

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
T.Trap = allTrapNumbersCol;
T.Strain = allStrainsCol;

writetable(T, tablefilename);
fprintf('Wrote %s with %d traps x %d timepoints.\n', tablefilename, size(T,1), numel(timeOffsets));
end
