function T = getCellAreasCummax(switchFrame, labels, tablefilename, varargin)

% getCellAreasCummax
% Cell area analysis using the maximum area of all segmented cells within each trap
%
% Input arguments:
% - switchFrame: scalar, frame index of media switch
% - labels: string array (2 elements) passed by wrapper
% - tablefilename: string/char, output CSV filename to write
%
% Name-value pairs:
% - 'TimeOffsets': numeric row vector, default [0 30 60 90]
%                  (offsets in frames or minutes, depending on your frame rate)
% - 'CodeDir': string/char, default 'code'
% - 'OutputDir': string/char, required by the wrapper
%
% Output:
%   T: table of max cell areas

%% Parse inputs

p = inputParser;
p.addParameter('TimeOffsets', [0 30 60 90]);
p.addParameter('CodeDir', 'code');
p.addParameter('OutputDir', '');
p.parse(varargin{:});

timeOffsets = p.Results.TimeOffsets;
outputDir = p.Results.OutputDir;

% Ensure row vector
timeOffsets = timeOffsets(:).';

%% Load experiment data from OutputDir

expInfoObj = loadExpInfo(outputDir);
posList = expInfoObj.getPositionList();
param = expInfoObj.getParameters();
nGrowthChannels = param.nGrowthChannels - length(param.emptyChannel);

%% Experimental settings

% Strain names
strainNames = [labels(1), labels(2)];

%% Main code

H = param.roiInPhase(4) + 1;
W = param.roiInPhase(3) + 1;
phaseChanName = expInfoObj.getChannelNames('phase');

% Determine how many frames (global across positions)
globalLastFrame = inf;
for pp = 1:length(posList)
    frameRange = expInfoObj.getRange(pp, phaseChanName);
    globalLastFrame = min(globalLastFrame, length(frameRange));
end

ycutoff = 150;

allCellAreas = cell(length(posList),1);
allStrains = cell(length(posList),1);
allAreaBeforeSwitch = cell(length(posList),1);

%% Define frame range 

% Maximum requested offset after switch
maxOffset = max(timeOffsets);

% Make sure don't request frames beyond what exists
maxOffsetAvailable = globalLastFrame - switchFrame;
if maxOffset > maxOffsetAvailable
    error('Requested TimeOffsets up to %d, but only %d frames are available after switchFrame.', ...
        maxOffset, maxOffsetAvailable);
end

% Number of frames before switch used to compute baseline 
preFrames = min(50, switchFrame - 1);

% Absolute frame indices used for analysis
frames = (switchFrame - preFrames):(switchFrame + maxOffset);

% Index of the switch frame within frames
idxSwitch = preFrames + 1;

%% Loop over positions

parfor pp = 1:length(posList)
    % Build trap mask
    trapLocs = expInfoObj.getChannelLocations(pp);
    trapLocs = trapLocs(2:end,:,1); 

    trapMask = zeros(H, W);
    for ii = 1:size(trapLocs,1)
        trapMask( ...
            trapLocs(ii,2):trapLocs(ii,2)+trapLocs(ii,4), ...
            trapLocs(ii,1):trapLocs(ii,1)+trapLocs(ii,3) ) = ii;
    end

    imDir = fullfile(outputDir, posList{pp}, "SegmentedPhase");
    imList = dir(fullfile(imDir, "*.tif*"));

    nFrames = numel(frames);
    A = zeros(nGrowthChannels, nFrames);
    nCells  = zeros(nGrowthChannels, nFrames);

    % Iterate over frames
    for i = 1:nFrames
        f = frames(i);

        % Skip if out of bounds
        if f < 1 || f > numel(imList)
            continue;
        end

        im = imread(fullfile(imDir, imList(f).name));

        props = regionprops(im, trapMask, "Area", "MaxIntensity", "BoundingBox");

        if isempty(props)

            % nothing segmented in this frame
            continue;
        end

        areas = vertcat(props.Area);
        ind = areas > 0;
        areas = areas(ind);
        traps = vertcat(props(ind).MaxIntensity);

        if f < switchFrame
            bb = vertcat(props(ind).BoundingBox);
            ind2 = bb(:,2) > ycutoff;
            areas = areas(ind2);
            traps = traps(ind2);
        end

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

    % Discard traps with fewer than 3 cells at the switch frame
    areaAfterSwitch(nCells(:, idxSwitch) < 3, :) = nan;

    % Select only the requested offsets
    offsetIdx = timeOffsets + 1;  % offset 0 -> col 1
    normalizedAreaAfterSwitch = areaAfterSwitch(:, offsetIdx);

    allCellAreas{pp} = normalizedAreaAfterSwitch;
    allAreaBeforeSwitch{pp} = areaBeforeSwitch;

    % first half positions -> labels(1), second half -> labels(2)
    if pp <= length(posList)/2
        strain = strainNames(1);
    else
        strain = strainNames(2);
    end
    allStrains{pp} = repmat(strain, size(normalizedAreaAfterSwitch,1), 1);
end

%% Collect across positions and build table

allCellAreas = cell2mat(allCellAreas);
allStrains = vertcat(allStrains{:});
allAreaBeforeSwitch = cell2mat(allAreaBeforeSwitch);

% Column names "0 min", "1 min", ..., "90 min"
timeVarNames = cellstr(strcat(string(timeOffsets), " min"));

% One column per time point
T = array2table(allCellAreas, 'VariableNames', timeVarNames);

% Add the extra columns
T.MeanAreaBeforeSwitch = allAreaBeforeSwitch;
T.Strain  = allStrains;

% Write with the filename
writetable(T, tablefilename);

end

