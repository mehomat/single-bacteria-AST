function [T] = getSingleLineageGrowthRatesInTBChip(outputDir, switchFrame, dt, strain, varargin)
% Computes lineage growth rates per chamber in giant chip (mother + descendants) and keeps only
% traps whose mother selection passes pre-switch growth QC
%
% Input:
% - outputDir: position
% - switchFrame: scalar, time point of switch to antibiotics
% - dt: scalar, time step
% - strain: string, name of the strain
%
% Name-value pairs:
% - 'Window': [NB NF] or scalar, default [10 0]
% - 'FitType': 'exp1' or 'poly1', default 'exp1'
% - 'Positions': numeric array, default []
%
% Output:
% - T: table including growth rates
%
% See also: getGrowthRateBeforeSwitch, getSingleLineageGrowthRate 
%% Parse inputs
p = inputParser;

p.addParameter('Window', [10 0], @(x) isnumeric(x) && (isscalar(x) || (numel(x)==2)));
p.addParameter('FitType', 'exp1', @(x) (ischar(x) || isstring(x)));
p.addParameter('Positions', [], @(x) isnumeric(x));

p.parse(varargin{:});

window = p.Results.Window;
fittype = p.Results.FitType;
posRange = p.Results.Positions;

%% Hard-coded parameters

minLength = 5; % min track length for growth rate fitting
maxGap = 5/dt; % max allowed gap in lineage tracking is 5 min
minCellGR = 0.002; % min allowed cell growth rate before adding AB
trackFrameRange = [switchFrame-10/dt, switchFrame+5/dt];
nChambers = 4;

%% Main code
allGrowthRates = [];
allFrames = [];
allTraps = [];
allPositions = [];
allCellIDs = [];

for chamberIdx=1:nChambers
    chamberOutputDir = fullfile(outputDir,['Chamber' num2str(chamberIdx)]);
    expInfoObj = loadExpInfo(chamberOutputDir);
    curPosRange = expInfoObj.getParameters('rangePositions');
    if isempty(posRange)
        validPosRange = curPosRange;
        validInd = 1:length(validPosRange);
    else
        [validPosRange,validInd] = intersect(curPosRange,posRange);
    end
    parfor pi = 1:length(validPosRange)
        pos = validPosRange(pi);
        posIdx = validInd(pi);

        % Load cell data for the current position 
        mCells = expInfoObj.getMCells(posIdx);

        % Compute single cell growth rate
        cellGR = nan(size(mCells));
        for i=1:length(mCells)
            if mCells(i).lifeTime > minLength
                cellGR(i) = getGrowthRateBeforeSwitch(mCells(i), switchFrame);
            end
        end

        isExplored = false(size(mCells));
        lastFrames = [mCells.lastFrame];
        trapGrowthRates = {};
        trapFrames = {};
        cellIDs = {};
        for i=1:length(mCells)
            c=mCells(i);
            if ~isExplored(i) && cellGR(i)>minCellGR
                [grs, frs] = getSingleLineageGrowthRate(c,window,dt,fittype);
        
                if isempty(frs)
                    continue
                end
        
                % Label explored cells
                linids = getCellLineage(c);
                idx = lastFrames(linids)<=frs(end);
                linids = linids(idx);
                isExplored(linids) = true;
                
                % Truncate at first big gap
                f = find(diff(frs) > maxGap, 1);
                if ~isempty(f)
                    frs = frs(1:f);
                    grs = grs(1:f);
                end
                 
                % Keep only traps covering full tracking window
                if ~isempty(frs) && frs(1) < trackFrameRange(1) && frs(end) > trackFrameRange(2)
                    trapFrames{end+1} = frs(:);
                    trapGrowthRates{end+1} = grs(:);
                    cellIDs{end+1} = c.id*ones(length(frs),1);
                end
                % end
            end
        end
        frames = cell2mat(trapFrames');
        growthrates = cell2mat(trapGrowthRates');
        cellIDs = cell2mat(cellIDs');
        allGrowthRates = [allGrowthRates; growthrates];
        allTraps = [allTraps; chamberIdx*ones(size(growthrates))];
        allCellIDs = [allCellIDs; cellIDs];
        allFrames = [allFrames; frames];
        allPositions = [allPositions; pos*ones(size(growthrates))];
    end
end


%% Build output table
allStrains = repmat(string(strain), size(allTraps));
allLineages = string(allPositions) + "_" + string(allTraps) + "_" +string(allCellIDs);
T = table(allGrowthRates, allFrames, allLineages,allStrains, ...
    'VariableNames', {'GrowthRate','Frame','Lineage','Strain'});