function T = getCellAreas(endFrame, label, tablefilename, posRange, outputDir, switchFrame, varargin)

% Extract cell area trajectories up to and including switchFrame.
%
% For each valid trap, follows the first mother cell and its descendants
% across all frames up to switchFrame, selecting at each division the
% daughter pushed deepest into the trap (away from the open end), based
% on growthDir.
%
% Areas are read directly from mCells.areas and stored unmodified
% (bad segmentations -> NaN)
%
% Input:
% - endFrame: scalar, frame index of media switch. Trajectories are
%   truncated to frames <= endFrame.
% - label: strain label(s). Either:
%     * scalar string/char (one strain) or
%     * 2-element string array / cell (two strains)
% - tablefilename: string/char, output CSV filename
% - posRange: positions to include. Either:
%     * numeric vector of position indices (one strain) or
%     * 1x2 cell array {positionsLabel1, positionsLabel2} (two strains)
%   The shape must match label.
% - outputDir: string/char, path to expInfoObj
% - switchFrame: scalar, frame index of media switch to get valid traps
%
% Name-value pairs:
% - dt: scalar, minutes per frame, default 1
%
% Output:
% - T: table

%% Parse inputs
p = inputParser;
p.addParameter('dt', 1, @(x) isnumeric(x) && isscalar(x));
p.parse(varargin{:});

dt = p.Results.dt;

%% Load experiment
expInfoObj = loadExpInfo(outputDir);
posList = expInfoObj.getPositionList();
param = expInfoObj.getParameters();
nGrowthChannels = param.nGrowthChannels - length(param.emptyChannel);
maxFrame = max(cellfun("length", expInfoObj.imRange));

% Handle empty channels
if isempty(param.emptyChannel)
    emptyChIdx = [];
else
    channelIndices = 1:param.nGrowthChannels;
    channelIndices(param.emptyChannel(2:end)) = [];
    emptyChIdx = find(channelIndices == param.emptyChannel(1));
end

% Normalize label / posRange into matching cell arrays
if iscell(posRange)
    posGroups = posRange;
else
    posGroups = {posRange};
end

labelList = string(label);
labelList = labelList(:).';

if numel(labelList) ~= numel(posGroups)
    error(['getCellAreasRaw: number of labels (%d) must match number ', ...
        'of position groups (%d).'], numel(labelList), numel(posGroups));
end

% Map each position index to its strain (empty if not assigned)
strainPerPos = strings(1, numel(posList));
for gi = 1:numel(posGroups)
    pp_g = unique(posGroups{gi}(:).');
    pp_g = pp_g(pp_g >= 1 & pp_g <= numel(posList));
    strainPerPos(pp_g) = labelList(gi);
end

%% Storage
allRows = {};

%% Loop over positions
for pp = 1:numel(posList)

    % Skip positions not assigned to any strain
    if strainPerPos(pp) == ""
        continue;
    end
    strain = strainPerPos(pp);

    % Get valid traps using the switch frame
    % require traps that were dividing before the switch
    S = getValidTraps(expInfoObj, pp, switchFrame, dt, ...
        'Traps', 1:nGrowthChannels, ...
        'RequireDividingBeforeSwitch', 1);

    if isempty(S)
        continue;
    end

    validTraps = [S([S.isValid]).trap]';
    if isempty(validTraps)
        continue;
    end

    % Load mCells for this position
    mCells = expInfoObj.getMCells(pp);
    if isempty(mCells)
        continue;
    end

    % Global trap ID offset
    trapIDoffset = (pp - 1) * nGrowthChannels;

    % Trap pixel areas for this position
    trapLocs = expInfoObj.getChannelLocations(pp);
    if ~isempty(emptyChIdx)
        trapLocs(emptyChIdx, :, :) = [];
    end
    trapLocs = trapLocs(:, :, 1);
    trapPixelArea = trapLocs(:, 3) .* trapLocs(:, 4); % width * height

    %% Loop over valid traps
    for ti = 1:numel(validTraps)

        trap = validTraps(ti);

        % Find S entry for this trap
        si = find([S.trap] == trap, 1);
        if isempty(si) || isempty(S(si).motherIds)
            continue;
        end

        % Growth direction for this trap
        growthDir = S(si).growthDir;

        % First mother cell in this trap
        motherIds = S(si).motherIds;
        firstMotherId = motherIds(1);
        firstMother = mCells(firstMotherId);

        %% Build trajectory: first mother + descendants chain
        allFrames = [];
        allAreas = [];
        allCellIds = [];

        c = firstMother;

        while ~isempty(c)

            nDet = numel(c.areas);
            cFrames = c.birthFrame : c.birthFrame + nDet - 1;
            cAreas = double(c.areas(:).');
            bad_seg = c.badSegmentations(:).';
            cAreas(bad_seg > 0) = NaN;

            allFrames = [allFrames, cFrames];
            allAreas = [allAreas, cAreas];
            allCellIds = [allCellIds, repmat(c.id, 1, nDet)];

            if isempty(c.descendants)
                break;
            elseif numel(c.descendants) == 1
                c = c.descendants(1);
            else
                dc1 = c.descendants(1);
                dc2 = c.descendants(2);

                idx1 = find(dc1.badSegmentations == 0, 1);
                idx2 = find(dc2.badSegmentations == 0, 1);

                y1 = NaN; if ~isempty(idx1), y1 = dc1.centroids(idx1, 2); end
                y2 = NaN; if ~isempty(idx2), y2 = dc2.centroids(idx2, 2); end

                if all(isnan([y1, y2]))
                    break;
                elseif growthDir < 0
                    % growthDir < 0 --> new cells push others toward smaller y,
                    % so the deepest (closed-end) cell has larger y
                    if y1 >= y2, c = dc1; else, c = dc2; end
                else
                    % growthDir > 0 --> new cells push others toward larger y,
                    % so the deepest (closed-end) cell has smaller y
                    if y1 <= y2, c = dc1; else, c = dc2; end
                end
            end

        end

        % Remove duplicate frames 
        [~, ia] = unique(allFrames, 'last');
        ia = sort(ia); % restore chronological order
        allFrames = allFrames(ia);
        allAreas = allAreas(ia);
        allCellIds = allCellIds(ia);

        % Truncate to frames <= endFrame
        keep = allFrames <= endFrame;
        allFrames = allFrames(keep);
        allAreas = allAreas(keep);
        allCellIds = allCellIds(keep);

        if isempty(allFrames)
            continue;
        end

        %% Cap area at trap pixel area (only non-NaN so it preserves NaN gaps)
        capArea = trapPixelArea(trap);
        nonNaN = ~isnan(allAreas);
        allAreas(nonNaN) = min(allAreas(nonNaN), capArea);

        %% Time vector
        time_min = (allFrames - 1) * dt;

        %% Store rows
        trapID = trapIDoffset + trap;

        for ri = 1:numel(allFrames)
            if isnan(allAreas(ri))
                continue;
            end

            row.TrapID = trapID;
            row.Position = pp;
            row.Trap = trap;
            row.Strain = strain;
            row.Frame = allFrames(ri);
            row.Time_min = time_min(ri);
            row.Area = allAreas(ri);
            row.TrapArea = trapPixelArea(trap);
            row.CellID = allCellIds(ri);
            allRows{end+1} = row;

        end

    end
end

%% Build table
if isempty(allRows)
    warning('getCellAreasRaw: no data collected.');
    T = table();
    writetable(T, tablefilename);
    return;
end

T = struct2table(vertcat(allRows{:}));
writetable(T, tablefilename);
fprintf('getCellAreasRaw: wrote %d rows to %s\n', height(T), tablefilename);

end