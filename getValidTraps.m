function S = getValidTraps(expInfoObj, pos, switchFrame, dt, varargin)

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
% - 'PreMarginMinutes': scalar minutes, default is 0
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
% S(ti).growthDir is the growth direction of the cells based on the sign of
% meanshifty

p = inputParser;
p.addParameter('Traps', [], @(x) isnumeric(x));
p.addParameter('YThresh', []);
p.addParameter('MinCellGR', 0.002, @(x) isnumeric(x) && isscalar(x));
p.addParameter('PreMarginMinutes', 0, @(x) isnumeric(x) && isscalar(x));
p.addParameter('TrackFrameRange', [switchFrame-10/dt, switchFrame+5/dt], @(x) isnumeric(x) && numel(x)==2);
p.addParameter('RequireDividingBeforeSwitch', false, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});

trapRange = p.Results.Traps;
yThresh = p.Results.YThresh;
minCellGR = p.Results.MinCellGR;
preMarginMinutes = p.Results.PreMarginMinutes;
trackFrameRange = p.Results.TrackFrameRange;
requireDividingBeforeSwitch = logical(p.Results.RequireDividingBeforeSwitch);

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
lifeTimes = [mCells.lifeTime];
isBadCellVals = [mCells.isBadCell];

% isBadCell filter
okBadCell = (isBadCellVals == 0 | isBadCellVals == 4);

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

% Save growth direction as an output
growthDir = sign(meanshifty);

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
S(nTraps) = struct('pos', [], 'trap', [], 'trapID', [], 'motherIds', [],  'isValid', [], 'growthDir', []);

trapIDs = (pos-1) * nGrowthChannels + trapRange;
for ti = 1:numel(trapRange)
    trap = trapRange(ti);

    % Anchor set, goodGrowth anchors the mother band
    indAnchor = (cellTraps == trap & goodGrowth & okBadCell);

    if any(indAnchor)
        if meanshifty > 0
            yCutOff = min(cellYcoord(indAnchor), [], 'omitnan') + dy;
            trapMotherCellIds = find(cellTraps == trap & cellYcoord < yCutOff);
            trapMotherCellIds = trapMotherCellIds(okBadCell(trapMotherCellIds));
        else
            yCutOff = max(cellYcoord(indAnchor), [], 'omitnan') - dy;
            trapMotherCellIds = find(cellTraps == trap & cellYcoord > yCutOff);
            trapMotherCellIds = trapMotherCellIds(okBadCell(trapMotherCellIds));
        end

        % Remove early born mothers that are not growing
        removeInd = birthFrames(trapMotherCellIds) < cutoffFrame & ~goodGrowth(trapMotherCellIds);
        trapMotherCellIds(removeInd) = [];

        % If mother cell has good enough GR before switch, but it has not
        % divided yet and one needs the last cell area before
        % division for baseline purposes
        if requireDividingBeforeSwitch

            trapMotherCellIds = trapMotherCellIds( ...
                isDividing(trapMotherCellIds) & (lastFrames(trapMotherCellIds) < switchFrame) );
        end

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
    S(ti).growthDir = growthDir; 
   
end

end
