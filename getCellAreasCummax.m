function T = getCellAreasCummax(switchFrame, labels, tablefilename, posRange, outputDir, varargin)

% getCellAreasCummax
% Cell area analysis using the maximum area of all segmented cells within each trap
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
%
% Output:
% - T: table of max cell areas

%% Parse inputs

p = inputParser;
p.addParameter('TimeOffsets', [0 30 60 90]);
p.parse(varargin{:});

timeOffsets = p.Results.TimeOffsets;

% Ensure row vector
timeOffsets = timeOffsets(:).';

%% Load experiment data from OutputDir

expInfoObj = loadExpInfo(outputDir);
posList = expInfoObj.getPositionList();
param = expInfoObj.getParameters();
nGrowthChannels = param.nGrowthChannels - length(param.emptyChannel);
if isempty(param.emptyChannel)
    emptyChIdx = [];
else
    channelIndices = 1:param.nGrowthChannels;
    channelIndices(param.emptyChannel(2:end)) = [];
    emptyChIdx = find(channelIndices == param.emptyChannel(1));
end

%% Experimental settings

% Strain names
strainNames = [labels(1), labels(2)];

% Position ranges
posLabel1 = unique(posRange{1}(:).');
posLabel2 = unique(posRange{2}(:).');

%% Main code

H = param.roiInPhase(4) + 1;
W = param.roiInPhase(3) + 1;

% Determine how many frames (global across positions)
globalLastFrame = min(cellfun("length",expInfoObj.imRange));

ycutoff = 150;

allCellAreas = cell(length(posList),1);
allStrains = cell(length(posList),1);
allAreaBeforeSwitch = cell(length(posList),1);
allTrapIDs = cell(length(posList), 1);

%% Define frame range 

% Maximum requested offset after switch
maxOffsetRequested = max(timeOffsets);

% How many frames are actually available after switchFrame
maxOffsetAvailable = globalLastFrame - switchFrame;

% If we don't have enough frames after the switch, just go to the last frame
maxOffset = min(maxOffsetRequested, maxOffsetAvailable);

% Trim offsets that exceed what's available
if maxOffset < maxOffsetRequested
    warning('Only %d frames available after switchFrame (requested up to %d). Trimming TimeOffsets.', ...
        maxOffsetAvailable, maxOffsetRequested);
    timeOffsets = timeOffsets(timeOffsets <= maxOffsetAvailable);
end

% Number of frames before switch used to compute baseline 
preFrames = min(50, switchFrame - 1);

% Absolute frame indices used for analysis (end at max available frame)
frames = (switchFrame - preFrames):(switchFrame + maxOffset);

% Index of the switch frame within frames
idxSwitch = preFrames + 1;


%% Loop over positions
parfor pp = 1:length(posList)

    % Get valid traps
    S = getValidTraps(expInfoObj, pp, switchFrame, 1, ...
        'Traps', 1:nGrowthChannels);

    validTraps = [S([S.isValid]).trap]';

    % Build trap mask
    trapLocs = expInfoObj.getChannelLocations(pp);
    if ~isempty(emptyChIdx)
        trapLocs(emptyChIdx,:,:)=[];
    end

    % Extract the x-coordinate
    trapLocs = trapLocs(:,:,1);
    %trapLocs = trapLocs(2:end,:,1); 

    trapMask = zeros(H, W);
    for ii = 1:size(trapLocs,1)
        trapMask( ...
            trapLocs(ii,2):trapLocs(ii,2)+trapLocs(ii,4), ...
            trapLocs(ii,1):trapLocs(ii,1)+trapLocs(ii,3) ) = ii;
    end

    imDir = fullfile(outputDir, posList{pp}, "SegmentedPhase");
    imList = dir(fullfile(imDir, "*.tif*"));

    nFrames = numel(frames);
    A = nan(nGrowthChannels, nFrames);
    nCells  = zeros(nGrowthChannels, nFrames);

    % Iterate over frames
    for i = 1:nFrames
        f = frames(i);

        % Skip if out of bounds
        if f < 1 || f > numel(imList)
            continue;
        end

        im = imread(fullfile(imDir, imList(f).name));

        % Mode based trap assignment
        props = regionprops(im, "Area", "PixelIdxList", "BoundingBox");

        if isempty(props)
            continue;
        end

        areasAll = vertcat(props.Area);

        % Compute trap id per object by majority vote of trapMask under the object
        trapsAll = nan(numel(props), 1);
        for k = 1:numel(props)
            vals = trapMask(props(k).PixelIdxList);
            vals = vals(vals > 0); % ignore outside-trap pixels
            if ~isempty(vals)
                trapsAll(k) = mode(vals); % most common trap id
            end
        end

        % Keep only valid objects (positive area AND assigned to a trap)
        keep = (areasAll > 0) & ~isnan(trapsAll);

        % Apply pre-switch ycutoff filter using the same keep indexing
        if f < switchFrame
            bbAll = vertcat(props.BoundingBox);
            keep = keep & (bbAll(:,2) > ycutoff);
        end

        areas = areasAll(keep);
        traps = trapsAll(keep);

        maxAreas = nan(nGrowthChannels, 1);
        for t = 1:nGrowthChannels
            indt = (traps == t);
            if any(indt)
                maxAreas(t) = max(areas(indt));
            end
            nCells(t,i) = sum(indt);
        end

        A(:,i) = maxAreas;
    end

    %% Compute pre-switch baseline and post-switch cummax

    % Frames before switch: columns 1:preFrames
    areaBeforeSwitch = round(mean(A(:, 1:preFrames), 2, 'omitnan'));

    % Frames from switch onwards
    A_afterSwitch   = A(:, idxSwitch:end);

    % Cumulative max over time after switch
    areaAfterSwitch = cummax(A_afterSwitch, 2);

    % Discard traps that are not valid by getValidTraps
    if isempty(validTraps)
        areaAfterSwitch(:,:) = nan;
    else
        invalidTraps = true(nGrowthChannels, 1);
        invalidTraps(validTraps) = false;
        areaAfterSwitch(invalidTraps, :) = nan;
    end

    % Select only the requested offsets
    offsetIdx = timeOffsets + 1;  % offset 0 -> col 1
    normalizedAreaAfterSwitch = areaAfterSwitch(:, offsetIdx);

    allCellAreas{pp} = normalizedAreaAfterSwitch;
    allAreaBeforeSwitch{pp} = areaBeforeSwitch;

    % Global Trap IDs for this position (pp)
    trapIDsGlobal = (pp-1)*nGrowthChannels + (1:nGrowthChannels).';
    allTrapIDs{pp} = trapIDsGlobal;

    % Assign strain by position ranges
    if ismember(pp, posLabel1)
        strain = strainNames(1); % labels(1)
    elseif ismember(pp, posLabel2)
        strain = strainNames(2); % labels(2)
    else
        strain = ""; % or unassigned
    end

    allStrains{pp} = repmat(strain, size(normalizedAreaAfterSwitch,1), 1);
    
end

%% Collect across positions and build table

allCellAreas = cell2mat(allCellAreas);
allStrains = vertcat(allStrains{:});
allAreaBeforeSwitch = cell2mat(allAreaBeforeSwitch);
allTrapIDs = vertcat(allTrapIDs{:});

% Column names 0 min, 1 min, ..., 90 min
timeVarNames = cellstr(strcat(string(timeOffsets), " min"));

% One column per time point
T = array2table(allCellAreas, 'VariableNames', timeVarNames);

% Add the extra columns
T.Trap = allTrapIDs;
T.MeanAreaBeforeSwitch = allAreaBeforeSwitch;
T.Strain = allStrains;

% Write with the filename
writetable(T, tablefilename);

end
