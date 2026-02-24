function T = getCellAreasCummaxUPD(switchFrame, labels, tablefilename, posRange, outputDir, varargin)

% Cell area analysis using max segmented object area in the mother region of each trap
% Baseline per trap = last mother cell's area on its last good frame BEFORE switchFrame
%
% Input arguments:
% - switchFrame: scalar, frame index of media switch
% - labels: string array (2 elements) 
% - tablefilename: string/char, output CSV filename to write
% - posRange: 1x2 cell array where:
%  posRange{1} = numeric vector of positions for labels(1)
%  posRange{2} = numeric vector of positions for labels(2)
% - outputDir: string/char, where the expInfoObj is loaded from
%
% Name-value pairs:
% - 'TimeOffsets': numeric row vector, default [0 30 60 90]
% - 'MotherFrac': fraction of the trap, from constriction side, that counts 
% as the mother region (along the y direction)
%
% Output:
% - T: table of max cell areas


%% Parse inputs
p = inputParser;

p.addParameter('TimeOffsets', [0 30 60 90]); % in frames

% MotherFrac defines region inside each trap along the y direction:
% - If constriction is at top of image coordinates --> mother region = top fraction of trap
% - If constriction is at bottom --> mother region = bottom fraction of trap
p.addParameter('MotherFrac', 1/3, @(x) isnumeric(x) && isscalar(x) && x>0 && x<=1);

p.parse(varargin{:});

timeOffsets = p.Results.TimeOffsets(:).'; % Force to be a row vector
motherFrac = p.Results.MotherFrac;

%% Load experiment data

expInfoObj = loadExpInfo(outputDir);
posList = expInfoObj.getPositionList();

% Parameters 
param = expInfoObj.getParameters();

% Number of traps per position excluding empty traps
nGrowthChannels = param.nGrowthChannels - length(param.emptyChannel);

% Handle empty channels by removing them from channelLocations
if isempty(param.emptyChannel)
    emptyChIdx = []; 
else
    channelIndices = 1:param.nGrowthChannels;
    channelIndices(param.emptyChannel(2:end)) = [];
    emptyChIdx = find(channelIndices == param.emptyChannel(1)); % the index in the "kept channel list"
end

%% Strain/position mapping

strainNames = [labels(1), labels(2)];

% Ensure position vectors are row vectors
posLabel1 = unique(posRange{1}(:).');
posLabel2 = unique(posRange{2}(:).');

%% Image geometry

% Coordinate conventions in images:
% - row index r corresponds to y-coordinate (downwards)
% - column index c corresponds to x-coordinate (to the right)
H = param.roiInPhase(4) + 1; % height in pixels (rows)
W = param.roiInPhase(3) + 1; % width in pixels (cols)

% Global last frame (minimum over positions because some positions might
% end earlier)
globalLastFrame = min(cellfun("length", expInfoObj.imRange));

%% Frame range

% Largest requested offset after switch
maxOffsetRequested = max(timeOffsets);

% Largest offset that still exists in the dataset
maxOffsetAvailable = globalLastFrame - switchFrame;

% Can only use what exists
maxOffset = min(maxOffsetRequested, maxOffsetAvailable);

% If asked for offsets beyond data availability, trim 
if maxOffset < maxOffsetRequested
    warning('Only %d frames available after switchFrame (requested up to %d). Trimming TimeOffsets.', ...
        maxOffsetAvailable, maxOffsetRequested);
    timeOffsets = timeOffsets(timeOffsets <= maxOffsetAvailable);
end

% Here frames(1) is switchFrame and frames(end) will be swf +
% maxOffset
frames = switchFrame:(switchFrame + maxOffset);

% Number of frames in analysis window
nFrames = numel(frames);

%% Preallocate outputs per position

% allCellAreas{pp} will contain a matrix size [nGrowthChannels x nRequestedOffsets]
allCellAreas = cell(length(posList),1); % (nTraps x nTimeOffsets)

% allStrains{pp} will contain a vector length nGrowthChannels of strain labels
allStrains = cell(length(posList),1);

% allAreaBeforeSwitch{pp} will contain baseline (one value per trap)
allAreaBeforeSwitch = cell(length(posList),1); % (nTraps x 1)

% allTrapIDs{pp} will contain global trap IDs for each trap in that position
allTrapIDs = cell(length(posList),1);

% Store trap inclusion counts per position
nNonEmptyTraps_perPos = nan(length(posList), 1);
nIncludedTraps_perPos = nan(length(posList), 1);

%% Loop over positions

parfor pp = 1:length(posList)

    % Included here to avoid warning from parfor when running
    growthDir = 0;

    % Get valid traps + mother IDs + growth direction
    S = getValidTraps(expInfoObj, pp, switchFrame, 1, ...
        'Traps', 1:nGrowthChannels, ...
        'RequireDividingBeforeSwitch', 1);

    if ~isempty(S)
        nNonEmptyTraps_perPos(pp) = S(1).nNonEmptyTraps;
        nIncludedTraps_perPos(pp) = S(1).nIncludedTraps;
    else
        nNonEmptyTraps_perPos(pp) = 0;
        nIncludedTraps_perPos(pp) = 0;
    end

    % Valid traps that passed QC
    validTraps = [S([S.isValid]).trap]';

    % growthDir is constant within a position
    if ~isempty(S)
        growthDir = S(1).growthDir; % +1 y increases over time, -1 y decreases over time
    end

    % Build trap mask + trap bounding boxes 
    % trapLocs is [nTraps x 4] with trapLocs(t,:) = [x0, y0, w, h]
    % x0,y0 = top left corner of the trap rectangle
    % w, h = width/height 
    trapLocs = expInfoObj.getChannelLocations(pp);

    % Remove empty channels from trapLocs so it aligns with 1:nGrowthChannels
    if ~isempty(emptyChIdx)
        trapLocs(emptyChIdx,:,:) = [];
    end

    % Keep only the geometry slice
    trapLocs = trapLocs(:,:,1);

    % trapMask is image sized matrix where each pixel value is
    % 0 if not inside any trap rectangle
    % t (1..nGrowthChannels) if inside trap t rectangle
    trapMask = zeros(H, W);

    for ii = 1:size(trapLocs,1)

        % trap rectangle geometry 
        x0 = trapLocs(ii,1); % column index of left edge
        y0 = trapLocs(ii,2); % row index of top edge
        ww = trapLocs(ii,3); % width
        hh = trapLocs(ii,4); % height

        % Fill rectangle with trap index
        trapMask(y0:y0+hh, x0:x0+ww) = ii;

    end

    % Load all mCells at this position
    mCells = expInfoObj.getMCells(pp);

    % Initialize baseline vector for traps
    areaBeforeSwitch = nan(nGrowthChannels, 1);

    % Loop over all traps
    for t = 1:nGrowthChannels

        % Find S entry for this trap
        si = find([S.trap] == t, 1, 'first');
        if isempty(si) || isempty(S(si).motherIds) % If no mother IDs
            continue;
        end

        % motherIds are indices into mCells
        mothers = S(si).motherIds(:);

        % Keep last mother before or at switch (already sorted by getValidTraps by birthFrame)
        motherId = mothers(end);

        % Compute baseline area from that mother track:
        % - Start from its last detection (last frame)
        % - If last detection is bad, step backwards until a good one
        % - Return corresponding mCell.areas value
        areaBeforeSwitch(t) = getLastGoodAreaBeforeOrAtFrame(mCells(motherId), switchFrame);

    end

    %% Iterate frames after switch and measure max area among objects in mother region

    imDir = fullfile(outputDir, posList{pp}, "SegmentedPhase");
    imList = dir(fullfile(imDir, "*.tif*"));

    % A(t,i) stores max area in trap t at analysis frame i 
    A = nan(nGrowthChannels, nFrames);

    for i = 1:nFrames

        % Absolute frame index in the experiment
        f = frames(i);

        % Ensure frame exists 
        if f < 1 || f > numel(imList)
            continue;
        end

        % Read segmented mask image for this frame
        im = imread(fullfile(imDir, imList(f).name));

        % Measure object areas and pixel index lists
        % props(k).Area is pixel area of object k
        % props(k).PixelIdxList are linear indices into the image matrix
        props = regionprops(im, "Area", "PixelIdxList");

        if isempty(props)
            continue;
        end

        % Areas for all objects
        areasAll = vertcat(props.Area);

        % trapsAll(k) = trap assignment of object k 
        trapsAll = nan(numel(props), 1);

        % yExtreme(k) = "closest to constriction" y pixel coordinate (row index) for object k
        % Depends on orientation
        % If constriction is at top --> closest pixel has min row index
        % If constriction is at bottom --> closest pixel has max row index
        yExtreme = nan(numel(props), 1);

        for k = 1:numel(props)

            % Linear indices of object pixels in the full image
            pix = props(k).PixelIdxList;

            % Map object pixels to trap ids
            vals = trapMask(pix);

            % Keep only pixels that lie inside any trap rectangle
            vals = vals(vals > 0);

            if isempty(vals)

                % Object does not overlap any trap region
                continue;
            end

            % Majority vote trap assignment
            % the trap id that occupies most pixels of the object
            trapsAll(k) = mode(vals);

            % Convert linear indices -> (row, col)
            % r are row indices (y in image coordinates)
            [r, ~] = ind2sub([H W], pix);

            if growthDir >= 0
                % If growthDir >= 0, treat constriction as at top 
                % closest to constriction pixel is then the min row
                yExtreme(k) = min(r);
            else
                % If growthDir < 0, treat constriction as at bottom 
                % closest to constriction pixel is then the max row
                yExtreme(k) = max(r);
            end
        end

        % Keep only objects that
        % have positive area
        % were assigned to a trap
        % have a defined yExtreme
        keep = (areasAll > 0) & ~isnan(trapsAll) & ~isnan(yExtreme);

        % Filtered vectors
        areas = areasAll(keep); % area per kept object
        traps = trapsAll(keep); % trap id per kept object
        yExt = yExtreme(keep); % closest to constriction y per kept object

        % Initiliaze store max area among objects in trap t that are within mother region
        maxAreas = nan(nGrowthChannels, 1);

        for t = 1:nGrowthChannels

            localIdx = find(traps == t);

            if isempty(localIdx)
                continue;
            end

            % Trap geometry:
            % trapLocs(t,:) = [x0, y0, w, h]
            % trapTop and trapBottom are in row coordinates
            y0 = trapLocs(t,2);
            hh = trapLocs(t,4);

            trapTop = y0; % smallest row index inside trap rectangle
            trapBottom = y0 + hh; % largest row index inside trap rectangle

            % Define mother region as a horizontal line within the trap
            if growthDir >= 0
                yThresh = trapTop + motherFrac * (trapBottom - trapTop);
                pass = yExt(localIdx) <= yThresh;
            else
                yThresh = trapBottom - motherFrac * (trapBottom - trapTop);
                pass = yExt(localIdx) >= yThresh;
            end

            % Apply mother region filter within this trap
            localIdxPass = localIdx(pass);

            % If anything passes, take the max area
            if ~isempty(localIdxPass)
                maxAreas(t) = max(areas(localIdxPass));
            end

        end

        % Store for this frame
        A(:, i) = maxAreas;
    end

    % cummax after switch 
    areaAfterSwitch = cummax(A, 2);

    % Discard traps not valid
    if isempty(validTraps)
        areaAfterSwitch(:,:) = nan;
    else
        invalidTraps = true(nGrowthChannels, 1);
        invalidTraps(validTraps) = false;
        areaAfterSwitch(invalidTraps, :) = nan;
        areaBeforeSwitch(invalidTraps) = nan;
    end

    % Select requested offsets
    offsetIdx = timeOffsets + 1;

    % Keep only requested time points
    normalizedAreaAfterSwitch = areaAfterSwitch(:, offsetIdx);

    % Store outputs for this position
    allCellAreas{pp} = normalizedAreaAfterSwitch;
    allAreaBeforeSwitch{pp} = areaBeforeSwitch;

    % Global trap IDs for this position
    trapIDsGlobal = (pp-1)*nGrowthChannels + (1:nGrowthChannels).';
    allTrapIDs{pp} = trapIDsGlobal;

    % Assign strain by position
    if ismember(pp, posLabel1)
        strain = strainNames(1);
    elseif ismember(pp, posLabel2)
        strain = strainNames(2);
    else
        strain = "";
    end

    % Replicate strain once per trap
    allStrains{pp} = repmat(strain, size(normalizedAreaAfterSwitch,1), 1);

end

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

%% Build + write summary CSV

Summary = table( ...
    strainNames(:), ...
    [totalNonEmpty_1; totalNonEmpty_2], ...
    [totalIncluded_1; totalIncluded_2], ...
    [pctIncluded_1; pctIncluded_2], ...
    'VariableNames', {'Strain','TotalNonEmptyTraps','TotalIncludedTraps','PercentIncludedTraps'});

summaryFilename = replace(string(tablefilename), ".csv", "_Summary.csv");
writetable(Summary, summaryFilename);

%% Collect and write

% Results over positions
allCellAreas = cell2mat(allCellAreas);
allStrains = vertcat(allStrains{:});
allAreaBeforeSwitch = cell2mat(allAreaBeforeSwitch);
allTrapIDs = vertcat(allTrapIDs{:});

timeVarNames = cellstr(strcat(string(timeOffsets), " min"));

% Table with one column per requested offset
T = array2table(allCellAreas, 'VariableNames', timeVarNames);

T.Trap = allTrapIDs;
T.LastMotherAreaBeforeSwitch = allAreaBeforeSwitch;
T.Strain = allStrains;

% Output CSV
writetable(T, tablefilename);

end

%% Helper to get last good area before a frame

function a = getLastGoodAreaBeforeOrAtFrame(mCell, frameCutoff)

% Returns the last good area for an mCell track on or before frameCutoff
% Good means badSegmentations == 0

a = NaN;

if ~isprop(mCell, 'areas') || isempty(mCell.areas)
    return;
end

nDet = numel(mCell.areas);

if isprop(mCell, 'badSegmentations') && ~isempty(mCell.badSegmentations)
    good = (mCell.badSegmentations == 0);
else
    good = true(1, nDet);
end

% Only keep detections that are at/before frameCutoff
lastDetByTime = min(nDet, frameCutoff - mCell.birthFrame + 1);
if lastDetByTime < 1
    return;
end

good = good(1:lastDetByTime);

if ~any(good)
    return;
end

lastIdx = find(good, 1, 'last');
a = double(mCell.areas(lastIdx));

end

