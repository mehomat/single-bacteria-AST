function [summary, perTrap, funnel] = getPreSwitchTrajectoryLoss(switchFrame, labels, posRange, outputDir, varargin)

% Account for every valid trap trajectory that does NOT end up as a clean,
% AB attributable trajectory in getCellAreasAlignedFillAreas and why
%
% Two orthogonal things are tracked per trap so nothing is collapsed away:
%
% 1. Trajectory fate is where the last valid (non NaN) area frame falls:
%
%       NoData: no valid area at all
%       LostBeforeABEffect: lastDataFrame < minFillFrame
%                           (never filled forward, a loss this early is
%                           not trusted as an AB response, so it is most
%                           likely a tracking/segmentation artifact)
%       Filled: minFillFrame <= lastDataFrame < maxFrame
%               (would be filled forward)
%       FullLength: lastDataFrame reaches the end of the experiment
%
% 2. Exclusion flags are reasons the main function drops a trap even when
%    the trajectory itself is complete:
%
%      FailsBaseline: no usable last-good mother area for baseline
%      FailsGenTime: growth rate fit failed (<15 good mother
%                    frames or non-positive slope) -> generationTime
%                    is NaN -> main function does continue
%
% A trap can be FullLength/Filled and STILL be excluded via FailsGenTime.
% That count is surfaced explicitly as nGenTimeFailButComplete.
%
% minFillFrame = switchFrame + round(MinMinutesPostSwitch / dt), matching
% getCellAreasAlignedFillAreas. The strictly-before-switch subset is kept
% as nEndsBeforeSwitch
%
% Trap selection and trajectory construction mirror
% getCellAreasAlignedFillAreas exactly, so counts are directly comparable
%
% Input:
% - switchFrame: scalar, frame index of media switch
% - labels: string array (2 elements)
% - posRange: 1x2 cell array {positionsLabel1, positionsLabel2}
% - outputDir: string/char, path to expInfoObj
%
% Name-value pairs:
% - dt: scalar, minutes per frame, default 1
% - MinMinutesPostSwitch: scalar, minutes after switch before fill-forward.
%                         Sets the LostBeforeABEffect cutoff, default 5
% - PreSwitchWindowMinutes: scalar, length (min) of the pre-switch window
%                         for the baseline loss rate. A lineage is counted
%                         (nValidLineages) only if it was present at the
%                         window start. That is, its first valid frame is 
%                         at or before switchFrame - this/dt and it had not ended
%                         before then. Of those, nLoss = lineages whose last
%                         valid frame falls inside [start, switchFrame).
%                         LossRate_pct = 100*nLoss/nValidLineages. Default 30.
% - SummaryFile: char/string path, if given, the summary table is
%                written there as CSV, default "" (no file)
% - PerTrapFile: char/string path, if given, the per-trap table is
%                written there as CSV, default "" (no file)
% - FunnelFile: char/string path, if given, the per-strain funnel
%               table is written there as CSV, default "" (no file)
% - ELN: char/string; if given, added as the first column
%        of all tables so per experiment files can be concatenated later, default ""
%
% Output:
% - summary: 1 row table of overall counts / percentages
% - perTrap: table, one row per valid trap, all flags + fate + reason
% - funnel: table, one row per strain (+ pooled ALL) giving the clean
%           top down waterfall: AllTraps -> NonEmpty -> Valid -> InAnalysis
%            -> KeptCompleteCurves, plus overlapping loss diagnostics

%% Parse inputs

p = inputParser;
p.addParameter('dt', 1, @(x) isnumeric(x) && isscalar(x));
p.addParameter('MinMinutesPostSwitch', 5, @(x) isnumeric(x) && isscalar(x));
p.addParameter('PreSwitchWindowMinutes', 30, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('SummaryFile', "", @(x) ischar(x) || isstring(x));
p.addParameter('PerTrapFile', "", @(x) ischar(x) || isstring(x));
p.addParameter('FunnelFile', "", @(x) ischar(x) || isstring(x));
p.addParameter('ELN', "", @(x) ischar(x) || isstring(x));
p.parse(varargin{:});

dt = p.Results.dt;
minMinutesPostSwitch = p.Results.MinMinutesPostSwitch;
preSwitchWindowMin = p.Results.PreSwitchWindowMinutes;
summaryFile = string(p.Results.SummaryFile);
perTrapFile = string(p.Results.PerTrapFile);
funnelFile = string(p.Results.FunnelFile);
elnTag = string(p.Results.ELN);

%% Load experiment
expInfoObj = loadExpInfo(outputDir);
posList = expInfoObj.getPositionList();
param = expInfoObj.getParameters();
nGrowthChannels = param.nGrowthChannels - length(param.emptyChannel);
maxFrame = max(cellfun("length", expInfoObj.imRange));

strainNames = [labels(1), labels(2)];
posLabel1 = unique(posRange{1}(:).');
posLabel2 = unique(posRange{2}(:).');

minFillFrame = switchFrame + round(minMinutesPostSwitch / dt);

% Pre-switch loss window: [windowStart, switchFrame). A trajectory is
% "lost in the window" if its last valid frame (the same lastDataFrame
% that would trigger fill forward in getCellAreasAlignedFillAreas) falls
% inside it
preSwitchWindowStart = max(1, switchFrame - round(preSwitchWindowMin / dt));

%% Storage
rows = {};
posStrainLog = strings(0, 1); % strain of each processed position
posAllLog = []; % growth channel traps per processed position
posNonEmptyLog = []; % non-empty traps per processed position

%% Loop over positions
for pp = 1:numel(posList)

    % Valid traps
    S = getValidTraps(expInfoObj, pp, switchFrame, dt, ...
        'Traps', 1:nGrowthChannels, ...
        'RequireDividingBeforeSwitch', 1);
    if isempty(S)
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

    % Funnel accounting that is recorded before skipping positions that 
    % have non empty traps but no valid ones
    posStrainLog(end+1, 1) = strain;
    posAllLog(end+1, 1) = nGrowthChannels;
    posNonEmptyLog(end+1, 1) = S(1).nNonEmptyTraps;

    validTraps = [S([S.isValid]).trap]';
    if isempty(validTraps)
        continue;
    end

    mCells = expInfoObj.getMCells(pp);
    if isempty(mCells)
        continue;
    end

    trapIDoffset = (pp - 1) * nGrowthChannels;

    %% Loop over valid traps
    for ti = 1:numel(validTraps)

        trap = validTraps(ti);

        si = find([S.trap] == trap, 1);
        if isempty(si) || isempty(S(si).motherIds)
            continue;
        end

        growthDir = S(si).growthDir;
        motherIds = S(si).motherIds;
        lastMotherId = motherIds(end);
        lastMother = mCells(lastMotherId);

        %% Baseline check (replicates main, on the ORIGINAL lastMother)
        good = (lastMother.badSegmentations == 0);
        if ~any(good)
            failsBaseline = true;
            baselineArea = NaN;
        else
            baselineArea = double(lastMother.areas(find(good, 1, 'last')));
            failsBaseline = isnan(baselineArea) || baselineArea <= 0;
        end

        %% Growth rate / generationTime check (replicates main)
        lmAreas = double(lastMother.areas(:).');
        lmBadSeg = lastMother.badSegmentations(:).';
        lmGood = (lmBadSeg == 0) & (lmAreas > 0) & ~isnan(lmAreas);
        nGoodMotherFrames = sum(lmGood);

        growthRate = NaN;
        generationTime = NaN;
        if nGoodMotherFrames >= 15
            lmTime = (find(lmGood) - 1) * dt;
            lmLnArea = log(lmAreas(lmGood));
            pFit = polyfit(lmTime(:), lmLnArea(:), 1);
            growthRate = pFit(1);
            if growthRate > 0
                generationTime = log(2) / growthRate;
            end
        end
        failsGenTime = isnan(generationTime);

        %% Build the full area trajectory across all mother cells
        allFrames = [];
        allAreas = [];

        for mi = 1:numel(motherIds)
            mc = mCells(motherIds(mi));

            nDet = numel(mc.areas);
            mcFrames = mc.birthFrame : mc.birthFrame + nDet - 1;
            mcAreas = double(mc.areas(:).');

            bad_seg = mc.badSegmentations(:).';
            mcAreas(bad_seg > 0) = NaN;

            allFrames = [allFrames, mcFrames];
            allAreas = [allAreas, mcAreas];
        end

        %% Follow the descendant chain from the lastMother
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
                allAreas = [allAreas,  cAreas];
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

        % Keep last occurrence
        [~, ia] = unique(allFrames, 'last');
        ia = sort(ia);
        allFrames = allFrames(ia);
        allAreas  = allAreas(ia);

        %% Trajectory fate from the last valid (non-NaN) frame
        validMask = ~isnan(allAreas);
        hasValidData = any(validMask);

        if ~hasValidData
            firstDataFrame = NaN;
            lastDataFrame = NaN;
            endsBeforeSwitch = false;
            endsBeforeABEffect = false;
            inWindow = false;
            lostInWindow = false;
            fate = "NoData";
        else
            firstDataFrame = min(allFrames(validMask));
            lastDataFrame = max(allFrames(validMask));
            endsBeforeSwitch = lastDataFrame < switchFrame;
            endsBeforeABEffect = lastDataFrame < minFillFrame;

            % Lineage is IN the window only if it was already present at the
            % window start. It started at or before switch-window, and had
            % not ended before the window start. Of those, "lost" = its last
            % valid frame (the fill trigger event) falls inside the window.
            inWindow = (firstDataFrame <= preSwitchWindowStart) && ...
                       (lastDataFrame  >= preSwitchWindowStart);
            lostInWindow = inWindow && (lastDataFrame < switchFrame);
            if endsBeforeABEffect
                fate = "LostBeforeABEffect";
            elseif lastDataFrame < maxFrame
                fate = "Filled";
            else
                fate = "FullLength";
            end
        end

        % Whether the main function would write any rows for this trap
        writesRows = hasValidData && ~failsBaseline && ~failsGenTime;

        % Primary exclusion reason, in the main function's check order
        if failsBaseline
            exclusionReason = "Baseline";
        elseif failsGenTime
            exclusionReason = "GenTime";
        elseif ~hasValidData
            exclusionReason = "NoTrajectoryData";
        else
            exclusionReason = "None";
        end

        %% Store
        r.TrapID = trapIDoffset + trap;
        r.Position = pp;
        r.Trap = trap;
        r.Strain = strain;
        r.FirstDataFrame = firstDataFrame;
        r.LastDataFrame = lastDataFrame;
        r.MinutesAfterSwitch = (lastDataFrame - switchFrame) * dt; % <0 => before switch
        r.nGoodMotherFrames = nGoodMotherFrames;
        r.GrowthRate = growthRate;
        r.GenerationTime_min = generationTime;
        r.BaselineArea = baselineArea;
        r.HasValidData = hasValidData;
        r.FailsBaseline = failsBaseline;
        r.FailsGenTime = failsGenTime;
        r.EndsBeforeSwitch = endsBeforeSwitch;
        r.EndsBeforeABEffect = endsBeforeABEffect;
        r.InWindow = inWindow;
        r.LostInWindow = lostInWindow;
        r.TrajectoryFate = fate;
        r.WritesRows = writesRows;
        r.ExclusionReason = exclusionReason;
        rows{end+1} = r;

    end
end

%% Assemble outputs
if isempty(rows)
    warning('getPreSwitchTrajectoryLoss: no valid traps found.');
    summary = table();
    perTrap = table();
    funnel = table();
    if summaryFile ~= "", writetable(summary, summaryFile); end
    if perTrapFile ~= "", writetable(perTrap, perTrapFile); end
    if funnelFile ~= "", writetable(funnel, funnelFile);  end
    return;
end

perTrap = struct2table(vertcat(rows{:}));

% Totals for the global summary, derived from the per position logs
nTrapsTotal = sum(posAllLog);
nNonEmptyTotal = sum(posNonEmptyLog);

nValid = height(perTrap);
fate = perTrap.TrajectoryFate;

nNoData = sum(fate == "NoData");
nLostBeforeABEffect = sum(fate == "LostBeforeABEffect");
nFilled = sum(fate == "Filled");
nFullLength = sum(fate == "FullLength");

nEndsBeforeSwitch = sum(perTrap.HasValidData & perTrap.EndsBeforeSwitch);
nFailsGenTime = sum(perTrap.FailsGenTime);
nFailsBaseline = sum(perTrap.FailsBaseline);
nGenTimeFailButComplete = sum(perTrap.FailsGenTime & ...
    (fate == "Filled" | fate == "FullLength"));
nInFinalAnalysis = sum(perTrap.WritesRows);

% Pre-switch window loss. Only lineages present at the
% window start are counted of which "lost" = ended inside the window
nValidLineages = sum(perTrap.InWindow);
nLoss = sum(perTrap.LostInWindow);

summary = struct2table(struct( ...
    'SwitchFrame', switchFrame, ...
    'PreSwitchWindowMin', preSwitchWindowMin, ...
    'nTrapsTotal', nTrapsTotal, ...
    'nNonEmptyTraps', nNonEmptyTotal, ...
    'nValidTraps', nValid, ...
    'pctValidOfNonEmpty', 100 * nValid / nNonEmptyTotal, ...
    'nValidLineages', nValidLineages, ...
    'nLoss', nLoss, ...
    'LossRate_pct', 100 * nLoss / nValidLineages, ...
    'nLostBeforeSwitch', nEndsBeforeSwitch, ...
    'pctLostBeforeSwitch', 100 * nEndsBeforeSwitch / nValid, ...
    'nLostBeforeABEffect', nLostBeforeABEffect, ...
    'pctLostBeforeABEffect', 100 * nLostBeforeABEffect / nValid, ...
    'nNoData', nNoData, ...
    'nFilled', nFilled, ...
    'nFullLength', nFullLength, ...
    'nFailsGenTime', nFailsGenTime, ...
    'pctFailsGenTime', 100 * nFailsGenTime / nValid, ...
    'nGenTimeFailButComplete', nGenTimeFailButComplete, ...
    'nFailsBaseline', nFailsBaseline, ...
    'nInFinalAnalysis', nInFinalAnalysis));

%% Per strain funnel (strict sequential waterfall, each trap removed once)

% Removals are applied IN ORDER to the valid traps:
%  1. ended before the switch frame  (incl. any with no valid data)
%  2. ended before switch+MinMinutesPostSwitch (only the EXTRA ones, i.e.
%     reached the switch but died within the window)
%  3. of those that reached switch+min, no generation-time fit
%  4. of the rest, no usable baseline
%
% What remains = complete trajectories kept in the analysis.
% NOTE "ended before switch" is a SUBSET of "ended before switch+min";
% steps 1 and 2 split that set so nothing is counted twice

strainList = unique(posStrainLog, 'stable');
strainList = [strainList; "ALL"];

hasV = perTrap.HasValidData;
ebs = perTrap.EndsBeforeSwitch;
eba = perTrap.EndsBeforeABEffect;
fgt = perTrap.FailsGenTime;
fbl = perTrap.FailsBaseline;

notReachSwitch = ~hasV | ebs; % step 1: never reached the switch frame
reachedSwitch = hasV & ~ebs;
beforeABextra = reachedSwitch & eba; % step 2: reached switch, died within window
complete = hasV & ~eba; % reached switch+min (= Filled/FullLength)

pct = @(a, b) 100 * a / b; % NaN for 0/0

funnelRows = {};
for k = 1:numel(strainList)
    s = strainList(k);
    if s == "ALL"
        pmask = true(size(posStrainLog));
        tmask = true(height(perTrap), 1);
    else
        pmask = (posStrainLog == s);
        tmask = (perTrap.Strain == s);
    end

    allTraps = sum(posAllLog(pmask));
    nonEmpty = sum(posNonEmptyLog(pmask));
    valid = sum(tmask);

    % Pre-switch window loss (headline), per strain
    nLin = sum(tmask & perTrap.InWindow);
    nLos = sum(tmask & perTrap.LostInWindow);

    remBeforeSwitch = sum(tmask & notReachSwitch);
    leftSwitch = valid - remBeforeSwitch;

    remBeforeAB = sum(tmask & beforeABextra);
    leftAB = leftSwitch - remBeforeAB; % == complete count

    remNoGenTime = sum(tmask & complete & fgt);
    leftGenTime = leftAB - remNoGenTime;

    remNoBaseline = sum(tmask & complete & ~fgt & fbl);
    kept = leftGenTime - remNoBaseline;

    f.Strain = s;
    f.AllTraps = allTraps;
    f.NonEmpty = nonEmpty;
    f.pctNonEmpty_ofAll = pct(nonEmpty, allTraps);
    f.Valid = valid;
    f.pctValid_ofNonEmpty = pct(valid, nonEmpty);
    f.nValidLineages = nLin;
    f.nLoss = nLos;
    f.LossRate_pct = pct(nLos, nLin);
    f.Rem_EndedBeforeSwitch = remBeforeSwitch;
    f.Left_ReachedSwitch = leftSwitch;
    f.pctLeft_ReachedSwitch = pct(leftSwitch, valid);
    f.Rem_BeforeABEffectExtra = remBeforeAB;
    f.Left_ReachedABWindow = leftAB;
    f.pctLeft_ReachedABWindow = pct(leftAB, valid);
    f.Rem_NoGenTimeFit = remNoGenTime;
    f.Left_AfterGenTime = leftGenTime;
    f.pctLeft_AfterGenTime = pct(leftGenTime, valid);
    f.Rem_NoBaseline = remNoBaseline;
    f.Kept = kept;
    f.pctKept_ofValid = pct(kept, valid);
    f.pctKept_ofNonEmpty = pct(kept, nonEmpty);
    funnelRows{end+1} = f;
end
funnel = struct2table(vertcat(funnelRows{:}));

%% Tag with ELN (so per-experiment files concatenate cleanly) and write CSVs
if elnTag ~= ""
    summary.ELN = elnTag;
    summary = movevars(summary, 'ELN', 'Before', 1);
    perTrap.ELN = repmat(elnTag, height(perTrap), 1);
    perTrap = movevars(perTrap, 'ELN', 'Before', 1);
    funnel.ELN = repmat(elnTag, height(funnel), 1);
    funnel = movevars(funnel, 'ELN', 'Before', 1);
end

if summaryFile ~= ""
    writetable(summary, summaryFile);
    fprintf('getPreSwitchTrajectoryLoss: wrote summary -> %s\n', summaryFile);
end
if perTrapFile ~= ""
    writetable(perTrap, perTrapFile);
    fprintf('getPreSwitchTrajectoryLoss: wrote per-trap -> %s\n', perTrapFile);
end
if funnelFile ~= ""
    writetable(funnel, funnelFile);
    fprintf('getPreSwitchTrajectoryLoss: wrote funnel -> %s\n', funnelFile);
end

%% Console report --> strict sequential waterfall, per strain and pooled (ALL)
fprintf('\n=== Trajectory loss funnel (switchFrame = %d, +%gmin window) ===\n', ...
    switchFrame, minMinutesPostSwitch);
for k = 1:height(funnel)
    fr = funnel(k, :);
    fprintf('\n[%s]\n', fr.Strain);
    fprintf('  All growth-channel traps:               %5d\n', fr.AllTraps);
    fprintf('  Non-empty (>=1 cell):                   %5d   (%.1f%% of all)\n', ...
        fr.NonEmpty, fr.pctNonEmpty_ofAll);
    fprintf('  Valid (passed getValidTraps):           %5d   (%.1f%% of non-empty)\n', ...
        fr.Valid, fr.pctValid_ofNonEmpty);
    fprintf('  >> %g-min pre-switch window:             %5d lineages, %d lost (%.1f%%)\n', ...
        preSwitchWindowMin, fr.nValidLineages, fr.nLoss, fr.LossRate_pct);
    fprintf('  start from valid:                       %5d   (100.0%% of valid)\n', fr.Valid);
    fprintf('   - ended before switch frame:         -%4d -> %5d   (%.1f%% of valid)\n', ...
        fr.Rem_EndedBeforeSwitch, fr.Left_ReachedSwitch, fr.pctLeft_ReachedSwitch);
    fprintf('   - ended before switch+%gmin (extra):  -%4d -> %5d   (%.1f%% of valid)\n', ...
        minMinutesPostSwitch, fr.Rem_BeforeABEffectExtra, fr.Left_ReachedABWindow, fr.pctLeft_ReachedABWindow);
    fprintf('   - no generation-time fit:            -%4d -> %5d   (%.1f%% of valid)\n', ...
        fr.Rem_NoGenTimeFit, fr.Left_AfterGenTime, fr.pctLeft_AfterGenTime);
    fprintf('   - no baseline:                       -%4d -> %5d   (%.1f%% of valid)\n', ...
        fr.Rem_NoBaseline, fr.Kept, fr.pctKept_ofValid);
    fprintf('  = KEPT in analysis:                     %5d   (%.1f%% of valid, %.1f%% of non-empty)\n', ...
        fr.Kept, fr.pctKept_ofValid, fr.pctKept_ofNonEmpty);
end
fprintf('\n');

end