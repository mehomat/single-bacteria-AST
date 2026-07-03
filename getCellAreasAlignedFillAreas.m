function T = getCellAreasAlignedFillAreas(switchFrame, labels, tablefilename, posRange, outputDir, varargin)

% Extract cell area trajectories aligned to the last mother cell
% division before AB introduction (t = 0).
%
% If the lineage divides again within 0.5 * generationTime after the
% switch frame, that division becomes the new t0. Chains until no
% further quick divisions occur within the window.
%
% Areas are read directly from mCells.areas
%
% Input:
% - switchFrame: scalar, frame index of media switch
% - labels: string array (2 elements)
% - tablefilename: string/char, output CSV filename
% - posRange: 1x2 cell array {positionsLabel1, positionsLabel2}
% - outputDir: string/char, path to expInfoObj
%
% Name-value pairs:
% - dt: scalar, minutes per frame, default 1
% - MinMinutesPostSwitch: scalar, minimum minutes after switch before
%   filling forward, default 5
% - CountDropThreshold: scalar in (0,1), fractional drop in cell count
%   from switchFrame to end of experiment required to classify as trap
%   emptied, default 0.5
% - CountCheckWindow: scalar, number of frames used to average counts
%   near end of experiment, default 10
% - CountSmoothWindow: scalar, number of preceding frames (inclusive)
%   used to smooth cell counts for the fill area calculation, default 5
%
% Output:
% - T: table

%% Parse inputs
p = inputParser;
p.addParameter('dt', 1, @(x) isnumeric(x) && isscalar(x));
p.addParameter('MinMinutesPostSwitch', 5, @(x) isnumeric(x) && isscalar(x));
p.addParameter('CountDropThreshold', 0.5, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('CountCheckWindow', 10, @(x) isnumeric(x) && isscalar(x));
p.addParameter('CountSmoothWindow', 5, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.parse(varargin{:});

dt = p.Results.dt;
minMinutesPostSwitch = p.Results.MinMinutesPostSwitch;
countDropThreshold = p.Results.CountDropThreshold;
countCheckWindow = p.Results.CountCheckWindow;
countSmoothWindow = p.Results.CountSmoothWindow;

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

strainNames = [labels(1), labels(2)];
posLabel1 = unique(posRange{1}(:).');
posLabel2 = unique(posRange{2}(:).');

%% Storage
allRows = {};

%% Loop over positions
for pp = 1:numel(posList)

    % Get valid traps
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

    % Strain label for this position
    if ismember(pp, posLabel1)
        strain = strainNames(1);
    elseif ismember(pp, posLabel2)
        strain = strainNames(2);
    else
        strain = "";
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

    % Build cell counts matrix for this position
    gc = [mCells.growthChannel];
    bf = [mCells.birthFrame];
    lf = [mCells.lastFrame];
    bad = [mCells.isBadCell];
    goodGlobal = (bad == 0) | (bad == 3) | (bad == 4);

    counts = zeros(nGrowthChannels, maxFrame, 'single');

    for ci = 1:numel(mCells)
        if ~goodGlobal(ci)
            continue;
        end

        trapIdx = gc(ci);
        seg = mCells(ci).badSegmentations;

        f1 = max(bf(ci), 1);
        f2 = min(lf(ci), maxFrame);
        if f2 < f1
            continue;
        end

        idx1 = f1 - bf(ci) + 1;
        idx2 = f2 - bf(ci) + 1;
        ok = (seg(idx1:idx2) == 0);
        frange = f1:f2;
        frange = frange(ok);

        counts(trapIdx, frange) = counts(trapIdx, frange) + 1;

    end

    % Smoothed counts (trailing mean) for fill area calculation
    smoothedCounts = zeros(size(counts), 'single');
    for ch = 1:nGrowthChannels
        for ff = 1:maxFrame
            f_start = max(1, ff - countSmoothWindow + 1);
            smoothedCounts(ch, ff) = mean(counts(ch, f_start:ff));
        end
    end

    %% Loop over valid traps
    for ti = 1:numel(validTraps)

        trap = validTraps(ti);

        % Find S entry for this trap
        si = find([S.trap] == trap, 1);
        if isempty(si) || isempty(S(si).motherIds)
            continue;
        end

        % Get growth direction for this trap from S
        growthDir = S(si).growthDir;

        % Mother IDs sorted by birthFrame
        motherIds = S(si).motherIds;

        % Last mother cell before switchFrame that divided
        lastMotherId = motherIds(end);
        lastMother = mCells(lastMotherId);

        % t = 0 frame = lastFrame of this mother cell (the division event)
        t0_frame = lastMother.lastFrame;

        % Baseline area = last good area of this mother at t0_frame
        good = (lastMother.badSegmentations == 0);
        if ~any(good)
            continue;
        end
        baselineArea = double(lastMother.areas(find(good, 1, 'last')));
        if isnan(baselineArea) || baselineArea <= 0
            continue;
        end

        %% Estimate growth rate from last mother cell before switch
        lmAreas = double(lastMother.areas(:).');
        lmBadSeg = lastMother.badSegmentations(:).';
        lmGood = (lmBadSeg == 0) & (lmAreas > 0) & ~isnan(lmAreas);

        growthRate = NaN;
        generationTime  = NaN;

        if sum(lmGood) >= 15 % Require at least 15 good segmentation frames for growth rate fit
            lmTime = (find(lmGood) - 1) * dt;
            lmLnArea = log(lmAreas(lmGood));

            pFit = polyfit(lmTime(:), lmLnArea(:), 1);
            growthRate = pFit(1);

            if growthRate > 0
                generationTime = log(2) / growthRate;
            end
        end

        % Skip trap if generation time could not be estimated
        if isnan(generationTime)
            continue;
        end

        %% Shift t0 forward if division occurs within 0.5 * genTime
        % after the SWITCH FRAME. Chains until no more quick
        % divisions within the window.
        halfGenFrames = round(0.5 * generationTime / dt);
        maxShiftFrame = switchFrame + halfGenFrames; % absolute deadline
        t0Cell = lastMother;

        shifting = true;
        while shifting
            shifting = false;

            if isempty(t0Cell.descendants)
                break;
            end

            % Pick the daughter following growth direction
            if numel(t0Cell.descendants) == 1
                daughter = t0Cell.descendants(1);
            else
                dc1 = t0Cell.descendants(1);
                dc2 = t0Cell.descendants(2);

                idx1 = find(dc1.badSegmentations == 0, 1);
                idx2 = find(dc2.badSegmentations == 0, 1);

                y1 = NaN; if ~isempty(idx1), y1 = dc1.centroids(idx1, 2); end
                y2 = NaN; if ~isempty(idx2), y2 = dc2.centroids(idx2, 2); end

                if all(isnan([y1, y2]))
                    break;
                elseif growthDir < 0
                    if y1 >= y2, daughter = dc1; else, daughter = dc2; end
                else
                    if y1 <= y2, daughter = dc1; else, daughter = dc2; end
                end
            end

            % Does this daughter divide (have descendants)?
            if isempty(daughter.descendants)
                break;
            end

            % Did the daughter divide within the allowed window?
            divisionFrame = daughter.lastFrame;
            if divisionFrame <= maxShiftFrame

                % Update t0 to this division
                t0_frame = divisionFrame;

                % Update baseline area
                dGood = (daughter.badSegmentations == 0);
                if any(dGood)
                    newBaseline = double(daughter.areas(find(dGood, 1, 'last')));
                    if ~isnan(newBaseline) && newBaseline > 0
                        baselineArea = newBaseline;
                    end
                end

                % Continue chaining from this daughter
                t0Cell = daughter;
                shifting = true;

            end
        end

        %% Build the full area trajectory across all mother cells
        allFrames = [];
        allAreas = [];
        allCellIds = [];

        for mi = 1:numel(motherIds)
            mc = mCells(motherIds(mi));

            nDet = numel(mc.areas);
            mcFrames = mc.birthFrame : mc.birthFrame + nDet - 1;
            mcAreas = double(mc.areas(:).');

            bad_seg = mc.badSegmentations(:).';
            mcAreas(bad_seg > 0) = NaN;

            allFrames = [allFrames, mcFrames];
            allAreas = [allAreas, mcAreas];
            allCellIds = [allCellIds, repmat(mc.id, 1, nDet)];
        end

        % Also include descendants after switchFrame - start from the
        % ORIGINAL lastMother (not t0Cell), because t0Cell may have been
        % shifted forward and still want all the data
        c = lastMother;
        visited = false;

        while ~isempty(c)

            if visited
                nDet = numel(c.areas);
                cFrames = c.birthFrame : c.birthFrame + nDet - 1;
                cAreas = double(c.areas(:).');
                bad_seg = c.badSegmentations(:).';
                cAreas(bad_seg > 0) = NaN;

                allFrames = [allFrames, cFrames];
                allAreas = [allAreas, cAreas];
                allCellIds = [allCellIds, repmat(c.id, 1, nDet)];
            end

            visited = true;

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
                    if y1 >= y2, c = dc1; else, c = dc2; end
                else
                    if y1 <= y2, c = dc1; else, c = dc2; end
                end
            end

        end

        % Remove duplicate frames (keep last occurrence = latest mother)
        [~, ia] = unique(allFrames, 'last');
        ia = sort(ia); % restore chronological order
        allFrames = allFrames(ia);
        allAreas = allAreas(ia);
        allCellIds = allCellIds(ia);

        if isempty(allFrames)
            continue;
        end

        % Last frame with a valid (non NaN) area
        validMask = ~isnan(allAreas);
        if ~any(validMask)
            continue;
        end
        lastDataFrame = max(allFrames(validMask));

        %% Fill forward + trap-emptied check
        minFillFrame = switchFrame + round(minMinutesPostSwitch / dt);
        isFilledForward = false;
        isTrapEmptied = false;
        fillFrames = [];

        if lastDataFrame < maxFrame && lastDataFrame >= minFillFrame

            % Trap level check that counts at switch vs mean count near end
            countAtSwitch = double(counts(trap, switchFrame));
            endFrames = max(1, maxFrame - countCheckWindow) : maxFrame;
            countEnd = mean(double(counts(trap, endFrames)), 'omitnan');

            if countAtSwitch > 0
                overallDropVal = (countAtSwitch - countEnd) / countAtSwitch;
                isTrapEmptied = overallDropVal >= countDropThreshold;
            else
                overallDropVal = NaN;
                isTrapEmptied = false;
            end

            % Fill forward with area based on smoothed cell count
            capArea = trapPixelArea(trap);
            fillFrames = (lastDataFrame + 1) : maxFrame;
            fillAreas = zeros(1, numel(fillFrames));

            for fi = 1:numel(fillFrames)
                nSmooth = smoothedCounts(trap, fillFrames(fi));
                if nSmooth >= 0.5
                    fillAreas(fi) = capArea / round(nSmooth);
                else
                    fillAreas(fi) = capArea; % trap empty -> full trap area
                end
            end

            allFrames = [allFrames,fillFrames];
            allAreas = [allAreas, fillAreas];
            allCellIds = [allCellIds, NaN(1, numel(fillFrames))];
            isFilledForward = true;

        end

        %% Cap area at trap pixel area (only non NaN, preserve NaN gaps)
        capArea = trapPixelArea(trap);
        nonNaN = ~isnan(allAreas);
        allAreas(nonNaN) = min(allAreas(nonNaN), capArea);

        %% Build raw nCells vector for all frames (real + filled)
        allNCells = zeros(1, numel(allFrames));
        for fi = 1:numel(allFrames)
            f = allFrames(fi);
            if f >= 1 && f <= maxFrame
                allNCells(fi) = double(counts(trap, f));
            end
        end

        %% Convert to time relative to t0_frame
        time_min = (allFrames - t0_frame) * dt;
        relArea = allAreas / baselineArea;
        lnRelArea = log(relArea);

        % Per-row filled flag (true only for the filled frames)
        nRealFrames = numel(time_min) - (isFilledForward * numel(fillFrames));
        isFilledRow = [false(1, nRealFrames), true(1, numel(time_min) - nRealFrames)];

        % Store rows
        trapID = trapIDoffset + trap;

        for ri = 1:numel(time_min)
            if isnan(lnRelArea(ri))
                continue;
            end

            row.TrapID  = trapID;
            row.Position = pp;
            row.Trap = trap;
            row.Strain = strain;
            row.Time_min = time_min(ri);
            row.RelArea = relArea(ri);
            row.LnRelArea = lnRelArea(ri);
            row.BaselineArea = baselineArea;
            row.TrapArea = trapPixelArea(trap);
            row.T0_frame = t0_frame;
            row.SwitchFrame = switchFrame;
            row.CellID = allCellIds(ri);
            row.NCells = allNCells(ri);
            row.GrowthRate = growthRate;
            row.GenerationTime_min = generationTime;
            row.IsFilledForward = isFilledForward;
            row.IsFilledRow = isFilledRow(ri);
            row.IsTrapEmptied = isTrapEmptied;
            allRows{end+1} = row;

        end

    end
end

%% Build table
if isempty(allRows)
    warning('getCellAreasAlignedFillAreas: no data collected.');
    T = table();
    writetable(T, tablefilename);
    return;
end

T = struct2table(vertcat(allRows{:}));
writetable(T, tablefilename);
fprintf('getCellAreasAlignedFillAreas: wrote %d rows to %s\n', height(T), tablefilename);

end