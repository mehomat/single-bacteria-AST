function S = getValidTraps(expInfoObj, pos, switchFrame, dt, varargin)
%
% Identify valid traps and return mother cell IDs per trap
% 
% Input:
% - expInfoObj: MATLAB object
% - pos: scalar, position
% - switchFrame: scalar, time point of switch to antibiotics
% - dt: scalar, time step
%
% Name-value pairs:
% - 'Traps': array of trap indices, default is all growth channels
% - 'YThresh': scalar or [], default is []
% - 'MinCellGR': scalar, default is 0.002
% - 'PreMarginMinutes': scalar minutes, default is 15
% - 'TrackFrameRange': [start end] frames, default is [switchFrame-10/dt, switchFrame+5/dt])
% - 'RequireDividing': logical, default is false
%
% Output:
% - S: structure array where
% S(ti).pos is the input position
% S(ti).trap is the trap index from trapRange(ti)
% S(ti).trapID is the global trap ID
% S(ti).motherIds is a column vector of indices into mCells for the trap
% mother cells found in that trap and sorted by birth frame
% S(ti).isValid is a logical flag and is only true if the trap passed the
% checks and if motherIds is not empty

p = inputParser;
p.addParameter('Traps', [], @(x) isnumeric(x));
p.addParameter('YThresh', []);
p.addParameter('MinCellGR', 0.002, @(x) isnumeric(x) && isscalar(x));
p.addParameter('PreMarginMinutes', 15, @(x) isnumeric(x) && isscalar(x));
p.addParameter('TrackFrameRange', [switchFrame-10/dt, switchFrame+5/dt], @(x) isnumeric(x) && numel(x)==2);
p.addParameter('RequireDividing', false, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});

trapRange = p.Results.Traps;
yThresh = p.Results.YThresh;
minCellGR = p.Results.MinCellGR;
preMarginMinutes = p.Results.PreMarginMinutes;
trackFrameRange = p.Results.TrackFrameRange;
requireDividing = logical(p.Results.RequireDividing);

param = expInfoObj.getParameters();
nGrowthChannels = param.nGrowthChannels - length(param.emptyChannel);
if isempty(trapRange)
    trapRange = 1:nGrowthChannels;
end

% Load cell data for this position
mCells = expInfoObj.getMCells(pos);
birthFrames = [mCells.birthFrame];
lastFrames = [mCells.lastFrame];
cellTraps = [mCells.growthChannel];

% Per cell properties
cellYcoord = nan(size(cellTraps));
cellLengths = nan(size(cellTraps));
shifty = zeros(size(cellTraps));
isDividing = false(size(cellTraps));
cellGR = nan(size(cellTraps));

for i = 1:numel(mCells)
    f = find(mCells(i).badSegmentations == 0, 1);
    if ~isempty(f)
        cellYcoord(i) = mCells(i).boundingBoxes(2, f);
        cellLengths(i) = mCells(i).boundingBoxes(4, f);
        shifty(i) = mean(diff(mCells(i).centroids(mCells(i).badSegmentations==0, 2)));
        if ~isempty(mCells(i).descendants)
            isDividing(i) = true;
        end
    end
    cellGR(i) = getGrowthRateBeforeSwitch(mCells(i), switchFrame);
end

% Orientation/geometry
dy = quantile(cellLengths(birthFrames < switchFrame), 0.1) / 2;
meanshifty = mean(shifty(lastFrames < switchFrame), 'omitnan');

if meanshifty < 0
    cellYcoord = cellYcoord + cellLengths;
end

% Optional y-threshold
if ~isempty(yThresh)
    if meanshifty < 0
        cellYcoord(cellYcoord > yThresh) = NaN;
    else
        cellYcoord(cellYcoord < yThresh) = NaN;
    end
end

% Growth QC values
goodGrowth = cellGR >= minCellGR;

% Pre-margin cutoff
preMarginFrames = round(preMarginMinutes / dt);
cutoffFrame = max(switchFrame - preMarginFrames, 1);

% Output
nTraps = numel(trapRange);
S(nTraps) = struct('pos', [], 'trap', [], 'trapID', [], 'motherIds', [],  'isValid', [] );

trapIDs = (pos-1) * nGrowthChannels + trapRange;
for ti = 1:numel(trapRange)
    trap = trapRange(ti);

    % Anchor set, goodGrowth anchors the mother band
    indAnchor = (cellTraps == trap & goodGrowth);
    if requireDividing
        indAnchor = indAnchor & isDividing;
    end

    if any(indAnchor)
        if meanshifty > 0
            yCutOff = min(cellYcoord(indAnchor), [], 'omitnan') + dy;
            trapMotherCellIds = find(cellTraps == trap & cellYcoord < yCutOff);
        else
            yCutOff = max(cellYcoord(indAnchor), [], 'omitnan') - dy;
            trapMotherCellIds = find(cellTraps == trap & cellYcoord > yCutOff);
        end

        % Remove early born mothers that are not growing
        removeInd = birthFrames(trapMotherCellIds) < cutoffFrame & ~goodGrowth(trapMotherCellIds);
        trapMotherCellIds(removeInd) = [];
    else
        trapMotherCellIds = [];
    end

    % Sort by birth frame
    trapBirthFrames = birthFrames(trapMotherCellIds);
    [trapBirthFrames, sortInd] = sort(trapBirthFrames);
    trapMotherCellIds = trapMotherCellIds(sortInd);

    if numel(trapBirthFrames) > numel(unique(trapBirthFrames))
        warning('Detected two trap mother cells in the frame in trap %d, pos %d', trap, pos)
    end

    % Require early start before tracking window
    isValid = true;
    if ~isempty(trapBirthFrames) && trapBirthFrames(1) > trackFrameRange(1)
        isValid = false;
        trapMotherCellIds = [];
    end

    % Store
    S(ti).pos = pos;
    S(ti).trap = trap;
    S(ti).trapID = trapIDs(ti);
    S(ti).motherIds = trapMotherCellIds(:); % column vector
    S(ti).isValid = isValid && ~isempty(trapMotherCellIds);
    
end

end
