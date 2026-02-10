function T = getTrapMotherCellLineageGrowthRatesSZ(expInfoObj, switchFrame, dt, strain, varargin)

% T = getTrapMotherCellLineageGrowthRates(expInfoObj,window,fittype,switchFrame,dt,strain,posRange)
% This function estimates the growth rate of the mother cells at the bottom of 
% mother-machine traps and their descendants for a strain treated with 
% antibiotics. Mother cells before the switch must be growing.
%
% Inputs:
% - expInfoObj: expInfo object with image analysis info
% - window: scalar or tuple
% - fittype: string/char, select from {'exp1','poly1'}
% - dt: scalar
%
% Name-value pairs:
% - 'Window': [NB NF] or scalar, default [10 0]
% - 'FitType': 'exp1' or 'poly1', default 'exp1'
% - 'Positions': scalar/array, default []
% - 'Traps': scalar/array, default []
% - 'YThresh': scalar or [], default []

%% Parse inputs
p = inputParser;

p.addParameter('Window', [10 0], @(x) isnumeric(x) && (isscalar(x) || (numel(x)==2)));
p.addParameter('FitType', 'exp1', @(x) (ischar(x) || isstring(x)));

p.addParameter('Positions', [], @(x) isnumeric(x));
p.addParameter('Traps', [], @(x) isnumeric(x));
p.addParameter('YThresh', []);
p.addParameter('Verbose', true);
p.parse(varargin{:});

window = p.Results.Window;
fittype = char(p.Results.FitType); % ensure char
posRange = p.Results.Positions;
trapRange = p.Results.Traps;
yThresh = p.Results.YThresh;

% Hard-coded parameters
minLength = 5; % min track length for growth rate fitting
maxGap = 5/dt; % max allowed gap in lineage tracking is 5 min
minCellGR = 0.002; % min allowed cell growth rate before adding AB
preMarginMinutes = 15; % How long before the switch one demands clean growth for mother cells
trackFrameRange = [switchFrame-10/dt, switchFrame+5/dt]; 

posList = expInfoObj.getPositionList();
if isempty(posRange), posRange = 1:length(posList); end
param = expInfoObj.getParameters();
nGrowthChannels = param.nGrowthChannels - length(param.emptyChannel);
if isempty(trapRange), trapRange = 1:nGrowthChannels; end

%% Main code

allGrowthRates = []; % growth rates
allFrames = []; % frames
allTraps = []; % traps
allStrains = [];

preMarginFrames = round(preMarginMinutes / dt);
cutoffFrame = max(switchFrame - preMarginFrames, 1);

for pi = 1:length(posRange)
    pos = posRange(pi);

    % Load cell data for the current position
    mCells = expInfoObj.getMCells(pos);
    birthFrames = [mCells.birthFrame];
    lastFrames = [mCells.lastFrame];
    cellTraps = [mCells.growthChannel];

    % Per-cell properties for mother selection
    cellYcoord = nan(size(cellTraps));
    cellLengths = nan(size(cellTraps));
    isDividing = false(size(cellTraps));
    shifty = zeros(size(cellTraps));
    cellGR = nan(size(cellTraps)); % growth rate before switch

    for i = 1:length(mCells)
        f = find(mCells(i).badSegmentations == 0, 1);
        if ~isempty(f)
            cellYcoord(i) = mCells(i).boundingBoxes(2, f);
            cellLengths(i) = mCells(i).boundingBoxes(4, f);
            shifty(i) = mean(diff(mCells(i).centroids(mCells(i).badSegmentations==0, 2)));

            if ~isempty(mCells(i).descendants)
                isDividing(i) = true;
            end
        end

        % Growth rate before the media switch for this cell
        cellGR(i) = getGrowthRateBeforeSwitch(mCells(i), switchFrame);

    end

    % Geometry/orientation
    dy = quantile(cellLengths(birthFrames < switchFrame), 0.1) / 2;
    meanshifty = mean(shifty(lastFrames < switchFrame), 'omitnan');
    % fprintf('%s(%d): meanshifty = %.2f\n', posList{pos}, pos, meanshifty);

    % If mother cells at bottom, shift y so bottom has larger y
    if meanshifty < 0
        cellYcoord = cellYcoord + cellLengths;
    end

    % Optional y-threshold discard (defined by the user) 
    if ~isempty(yThresh)
        if meanshifty < 0
            cellYcoord(cellYcoord > yThresh) = NaN;
        else
            cellYcoord(cellYcoord < yThresh) = NaN;
        end
    end

    % Logical mask for "good growth before switch"
    goodGrowth = cellGR >= minCellGR;

    posFrames = [];
    posGrowthRates = [];
    posTraps = [];
    parfor ti = 1:length(trapRange)
        trap = trapRange(ti);

        % Find trap mother cells and sort them by the birth frame
        % Find trap mother cells that grow before media switch
        indGoodGRTrap = (cellTraps == trap & goodGrowth);  % Same logic as old code

        if any(indGoodGRTrap)
            if meanshifty > 0
                % Mother cells at the top
                yCutOff = min(cellYcoord(indGoodGRTrap), [], 'omitnan') + dy; % Identifies "mother band"
                trapMotherCellIds = find(cellTraps == trap & cellYcoord < yCutOff);
            else
                % Mother cells at the bottom
                yCutOff = max(cellYcoord(indGoodGRTrap), [], 'omitnan') - dy;
                trapMotherCellIds = find(cellTraps == trap & cellYcoord > yCutOff);
            end
            % Discard trap mother cells that do not grow before the switch
            removeInd = birthFrames(trapMotherCellIds)<cutoffFrame & ~goodGrowth(trapMotherCellIds);
            trapMotherCellIds(removeInd)=[];
        else
            trapMotherCellIds = [];
        end

        % Sort by birth frame
        [trapBirthFrames, sortInd] = sort(birthFrames(trapMotherCellIds));

        if numel(trapBirthFrames) > numel(unique(trapBirthFrames))
            warning('Detected two trap mother cells in the frame in trap %d, pos %d', trap, pos)
        end

        trapMotherCellIds = trapMotherCellIds(sortInd);


        % Require that the earliest mother in this trap starts before the tracking window
        if ~isempty(trapBirthFrames) && trapBirthFrames(1) > trackFrameRange(1)
            trapMotherCellIds = [];
        end

        % Lineage growth rate computation
        trapGrowthRates = [];
        trapFrames = [];
        for i = 1:length(trapMotherCellIds)
            cid = trapMotherCellIds(i);
            c = mCells(cid);

            [areas, badSegmentations] = getSumOfLineage(c);
            frs = c.birthFrame : c.birthFrame + length(areas) - 1;

            if ~isempty(trapFrames)
                firstFrameIndex = find(frs > trapFrames(end) - window(1), 1);
                areas = areas(firstFrameIndex:end);
                badSegmentations = badSegmentations(firstFrameIndex:end);
                frs = frs(firstFrameIndex:end);
            end

            if length(frs) > minLength
                t = frs * dt;
                grs = movgrowthrate2(t, areas, badSegmentations, window, fittype);

                % Skip overlapping frames
                if ~isempty(trapFrames)
                    firstFrameIndex = find(frs > trapFrames(end), 1);
                    grs = grs(firstFrameIndex:end);
                    frs = frs(firstFrameIndex:end);
                end

                trapGrowthRates = [trapGrowthRates; grs];
                trapFrames = [trapFrames; frs(:)];
            end
        end

        % Clean up and truncate long gaps
        indFinite = isfinite(trapGrowthRates);
        trapGrowthRates = trapGrowthRates(indFinite);
        trapFrames = trapFrames(indFinite);

        f = find(diff(trapFrames) > maxGap, 1);
        if ~isempty(f)
            trapFrames = trapFrames(1:f);
            trapGrowthRates = trapGrowthRates(1:f);
        end

        % Keep only traps covering full tracking window
        if ~isempty(trapFrames) && trapFrames(1) < trackFrameRange(1) && trapFrames(end) > trackFrameRange(2)

            % DEBUG (Print mother cell IDs used for this trap)
            %idStr = sprintf('%d ', trapMotherCellIds);
            %fprintf('Accepted trap: pos %d, trap %d (trapID %d). Mother cell IDs: [%s]\n', ...
                    %pos, trap, (pos-1)*nGrowthChannels + trap, strtrim(idStr));

            % Store results
            trapID = (pos-1) * nGrowthChannels + trap;
            posGrowthRates = [posGrowthRates; trapGrowthRates];
            posFrames = [posFrames; trapFrames];
            posTraps = [posTraps; repmat(trapID, size(trapGrowthRates))];
        end
    end

    allGrowthRates = [allGrowthRates; posGrowthRates];
    allFrames = [allFrames; posFrames];
    allTraps = [allTraps; posTraps];
    
end

if ~isempty(strain)
    allStrains = repmat(strain, size(allTraps));
    T = table(allGrowthRates, allFrames, allTraps, allStrains, ...
        'VariableNames', {'GrowthRate','Frame','Trap','Strain'});
else
    T = table(allGrowthRates, allFrames, allTraps, ...
        'VariableNames', {'GrowthRate','Frame','Trap'});
end
end
