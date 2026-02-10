function T = getTrapMotherCellLineageGrowthRatesVaryingWindow(expInfoObj, switchFrame, dt, strain, varargin)

% Computes lineage growth rates per trap (mother + descendants), using:
% - WindowBefore for frames <= switchFrame
% - WindowAfterSwitch for frames > switchFrame
%
% Input:
% - expInfoObj: MATLAB object
% - switchFrame: scalar, time point of switch to antibiotics
% - dt: scalar, time step
% - strain: string, name of the strain
%
% Name-value pairs:
% - 'WindowBefore': [NB NF] or scalar, default [10 0]
% - 'WindowAfterSwitch': [NB NF] or scalar, default [10 0]
% - 'FitType': 'exp1' or 'poly1', default 'exp1'
% - 'Positions': numeric array, default []
% - 'Traps': numeric array, default []
% - 'YThresh': scalar or [], default []
%
% Output:
% - T: table including growth rates

%% Parse inputs
p = inputParser;

p.addParameter('WindowBefore', [10 0], @(x) isnumeric(x) && (isscalar(x) || numel(x)==2));
p.addParameter('WindowAfterSwitch', [10 0], @(x) isnumeric(x) && (isscalar(x) || numel(x)==2));
p.addParameter('FitType', 'exp1', @(x) (ischar(x) || isstring(x)));
p.addParameter('Positions', [], @(x) isnumeric(x));
p.addParameter('Traps', [], @(x) isnumeric(x));
p.addParameter('YThresh', []);

p.parse(varargin{:});

Kpre  = p.Results.WindowBefore;
Kpost = p.Results.WindowAfterSwitch;
fittype = char(p.Results.FitType);
posRange = p.Results.Positions;
trapRange = p.Results.Traps;
yThresh = p.Results.YThresh;

%% Hard-coded parameters
minLength = 5; % min track length for growth rate fitting
maxGap = 5/dt; % max allowed gap in lineage tracking is 5 min
minCellGR = 0.002; % min allowed cell growth rate before adding AB
preMarginMinutes = 15; % how long before switch one demands clean growth
trackFrameRange = [switchFrame-10/dt, switchFrame+5/dt];

posList = expInfoObj.getPositionList();
if isempty(posRange), posRange = 1:length(posList); end

param = expInfoObj.getParameters();
nGrowthChannels = param.nGrowthChannels - length(param.emptyChannel);
if isempty(trapRange), trapRange = 1:nGrowthChannels; end

%% Main code
allGrowthRates = [];
allFrames = [];
allTraps = [];

for pi = 1:length(posRange)
    pos = posRange(pi);

    % Load cell data for the current position
    mCells = expInfoObj.getMCells(pos);

    % Valid traps + mother IDs
    S = getValidTraps(expInfoObj, pos, switchFrame, dt, ...
        'Traps', trapRange, ...
        'YThresh', yThresh, ...
        'MinCellGR', minCellGR, ...
        'PreMarginMinutes', preMarginMinutes, ...
        'TrackFrameRange', trackFrameRange, ...
        'RequireDividing', false);

    validMask = [S.isValid];
    motherIdCells = {S.motherIds};
    trapIDs = [S.trapID];
    trapNums = [S.trap];

    % Per position results
    posFrames = [];
    posGrowthRates = [];
    posTraps = [];

    parfor ti = 1:numel(trapNums)
        if ~validMask(ti)
            continue;
        end

        trapID = trapIDs(ti);
        trapMotherCellIds = motherIdCells{ti};

        trapGrowthRates = [];
        trapFrames = [];

        for ii = 1:numel(trapMotherCellIds)
            cid = trapMotherCellIds(ii);
            c = mCells(cid);

            [areas, badSegmentations] = getSumOfLineage(c);
            frs = c.birthFrame : (c.birthFrame + length(areas) - 1);

            % Avoid too much overlap with previous segment
            if ~isempty(trapFrames)
                firstFrameIndex = find(frs > trapFrames(end) - Kpre(1), 1);
                if isempty(firstFrameIndex)
                    continue;
                end
                areas = areas(firstFrameIndex:end);
                badSegmentations = badSegmentations(firstFrameIndex:end);
                frs = frs(firstFrameIndex:end);
            end

            if numel(frs) > minLength
                t = frs * dt;

                % Compute pre/post GRs and stitch
                gr_pre = movgrowthrate2(t, areas, badSegmentations, Kpre,  fittype);
                gr_post = movgrowthrate2(t, areas, badSegmentations, Kpost, fittype);

                grs = stitchPrePost(gr_pre, gr_post, frs, switchFrame);

                % Skip overlap frames entirely
                if ~isempty(trapFrames)
                    firstFrameIndex = find(frs > trapFrames(end), 1);
                    if isempty(firstFrameIndex)
                        continue;
                    end
                    grs = grs(firstFrameIndex:end);
                    frs = frs(firstFrameIndex:end);
                end

                trapGrowthRates = [trapGrowthRates; grs(:)];
                trapFrames = [trapFrames; frs(:)];
            end
        end

        % Clean up
        indFinite = isfinite(trapGrowthRates);
        trapGrowthRates = trapGrowthRates(indFinite);
        trapFrames = trapFrames(indFinite);

        % Truncate at first big gap
        f = find(diff(trapFrames) > maxGap, 1);
        if ~isempty(f)
            trapFrames = trapFrames(1:f);
            trapGrowthRates = trapGrowthRates(1:f);
        end

        % Keep only traps covering full tracking window
        if ~isempty(trapFrames) && trapFrames(1) < trackFrameRange(1) && trapFrames(end) > trackFrameRange(2)
            posGrowthRates = [posGrowthRates; trapGrowthRates];
            posFrames = [posFrames; trapFrames];
            posTraps = [posTraps; repmat(trapID, size(trapGrowthRates))];
        end
    end

    allGrowthRates = [allGrowthRates; posGrowthRates];
    allFrames = [allFrames; posFrames];
    allTraps = [allTraps; posTraps];
end

% Build output table
if ~isempty(strain)
    allStrains = repmat(string(strain), size(allTraps));
    T = table(allGrowthRates, allFrames, allTraps, allStrains, ...
        'VariableNames', {'GrowthRate','Frame','Trap','Strain'});
else
    T = table(allGrowthRates, allFrames, allTraps, ...
        'VariableNames', {'GrowthRate','Frame','Trap'});
end

end

function gr = stitchPrePost(gr_pre, gr_post, frs, switchFrame)

% Force column vectors 
gr_pre = gr_pre(:);
gr_post = gr_post(:);
frs = frs(:);

% Hard cut (pre for frames <= switchFrame, post for frames > switchFrame)
gr = gr_pre;
idx = frs > switchFrame;
gr(idx) = gr_post(idx);

end


