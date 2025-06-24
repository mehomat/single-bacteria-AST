function T = getTrapMotherCellLineageGrowthRates(expInfoObj,window,fittype,switchFrame,dt,strain,posRange)
% T = getTrapMotherCellLineageGrowthRates(expInfoObj,window,fittype,switchFrame,dt,strain,posRange)
% This function estimates the growth rate of the mother cells at the bottom of 
% mother-machine traps and their descendants for a strain treated with 
% antibiotics
%
% Input:
%   expInfoObj - expInfo object
%   window - if scalar K, the window is centered and has length K frames.
%            When K is even, the window is centered about the current and 
%            previous elements of X.
%            If tuple (NB,NF), the window uses previous NB frames, the 
%            current frame, and the next NF frames, i.e., in total NB+NF+1.
%   fittype - either 'fit1' (exponential fit) or 'poly1' (linear fit)
%   switchFrame - frame index, where the first media switch occurs
%   dt - frame rate, in minutes.
%   strain - string with the name of the string. Optional.
%   posRange - range of positions with the strain. Optional, use all
%              positions by default.
%
% Output: table T with columns 'Frame', 'Trap', 'GrowthRate', 'Strain'.
% The trap ID for trap t in position p is (p-1)*nTraps+t.

%parse parameters
if nargin<6, strain = "";end
if nargin<7, posRange = [];end

%hard-coded parameters
minLength = 5; % min track length for growth rate fitting
maxGap = 5/dt; %  max allowed gap in cell lineage tracking is 5 min
% track cells for more than 10 min before adding antibiotics and 
% more than 5 min after adding antibiotics
trackFrameRange = [switchFrame-10/dt, switchFrame+5/dt]; 

posList = expInfoObj.getPositionList();
if isempty(posRange), posRange=1:length(posList);end
param = expInfoObj.getParameters();
nGrowthChannels = param.nGrowthChannels-length(param.emptyChannel);


allGrowthRates = []; % growth rates
allFrames = []; % frames
allTraps = []; % traps
allStrains = [];
for pi=1:length(posRange)
    pos = posRange(pi);

    % Load cell data for the current position
    mCells = expInfoObj.getMCells(pos);
    birthFrames = [mCells.birthFrame];
    lastFrames = [mCells.lastFrame];
    cellTraps = [mCells.growthChannel];


    % Identify trap mother cell locations
    cellYcoord = nan(size(cellTraps));
    cellLengths = nan(size(cellTraps));
    isDividing = false(size(cellTraps));
    shifty = zeros(size(cellTraps));
    for i=1:length(mCells)
        f = find(mCells(i).badSegmentations==0,1);
        if ~isempty(f)
            cellYcoord(i)=mCells(i).boundingBoxes(2,f);
            cellLengths(i)=mCells(i).boundingBoxes(4,f);
            shifty(i)=mean(diff(mCells(i).centroids(mCells(i).badSegmentations==0,2)));
            if ~isempty(mCells(i).descendants)
                isDividing(i)=true;
            end
        end
    end
    dy=quantile(cellLengths(birthFrames<switchFrame),0.1)/2;
    meanshifty = mean(shifty(lastFrames<switchFrame),'omitnan');
    fprintf('%s(%d): meanshifty = %.2f\n',posList{pos},pos,meanshifty);
    if meanshifty<0 % trap mother cells located at the bottom of images
        cellYcoord = cellYcoord+cellLengths;
    end
    
    %
    posFrames = [];
    posGrowthRates = [];
    posTraps = [];
    parfor trap = 1:nGrowthChannels
    
        % find trap mother cells and sort them by the birth frame
        ind = cellTraps==trap & isDividing;
        if any(ind)
            if meanshifty>0
                yCutOff = min(cellYcoord(cellTraps==trap & isDividing))+dy;
                trapMotherCellIds = find(cellTraps==trap & cellYcoord<yCutOff);
            else
                yCutOff = max(cellYcoord(cellTraps==trap & isDividing))-dy;
                trapMotherCellIds = find(cellTraps==trap & cellYcoord>yCutOff);
            end
        else
            trapMotherCellIds = [];
        end

        [trapBirthFrames,sortInd] = sort(birthFrames(trapMotherCellIds));
        if numel(trapBirthFrames)>numel(unique(trapBirthFrames))
             warning('Detected two trap mother cells in the same frame in trap %d, pos %d',trap,pos)
        end
        trapMotherCellIds = trapMotherCellIds(sortInd);
        if ~isempty(trapBirthFrames) && trapBirthFrames(1)>trackFrameRange(1)
            trapMotherCellIds = [];
        end
        
        trapGrowthRates = [];
        trapFrames = [];
        for i=1:length(trapMotherCellIds)
            cid = trapMotherCellIds(i);
            c = mCells(cid);
            [areas, badSegmentations] = getSumOfLineage(c);
            frs = c.birthFrame : c.birthFrame+length(areas)-1;
            if ~isempty(trapFrames)
                firstFrameIndex = find(frs>trapFrames(end)-window(1),1);
                areas = areas(firstFrameIndex:end);
                badSegmentations = badSegmentations(firstFrameIndex:end);
                frs = frs(firstFrameIndex:end);
            end
            if length(frs)>minLength
                t = frs*dt;
                grs = movgrowthrate2(t,areas,badSegmentations,window,fittype);
                % skip overlap frames
                if ~isempty(trapFrames)
                    firstFrameIndex = find(frs>trapFrames(end),1);
                    grs = grs(firstFrameIndex:end);
                    frs = frs(firstFrameIndex:end);
                end
                trapGrowthRates = [trapGrowthRates; grs];
                trapFrames  = [trapFrames; frs'];
            end
        end

        ind = isfinite(trapGrowthRates);
        trapGrowthRates = trapGrowthRates(ind);
        trapFrames = trapFrames(ind);
        % if a cell lineage has a large gap, keep the first part
        f = find(diff(trapFrames)>maxGap,1);
        if ~isempty(f)
            trapFrames = trapFrames(1:f);
            trapGrowthRates = trapGrowthRates(1:f);
        end
        if ~isempty(trapFrames) && trapFrames(1)<trackFrameRange(1) && trapFrames(end)>trackFrameRange(2)
            posGrowthRates = [posGrowthRates; trapGrowthRates];
            posFrames = [posFrames; trapFrames];
            trapID = (pos-1)*nGrowthChannels+trap;
            posTraps = [posTraps; repmat(trapID,size(trapGrowthRates))];
        end
    end
    allGrowthRates = [allGrowthRates; posGrowthRates];
    allFrames = [allFrames; posFrames];
    allTraps = [allTraps; posTraps];
end

if ~isempty(strain)
    allStrains = [allStrains; repmat(strain,size(allTraps))];
    T = table(allGrowthRates,allFrames,allTraps,allStrains,...
        'VariableNames',{'GrowthRate','Frame','Trap','Strain'});
else
    T = table(allGrowthRates,allFrames,allTraps,...
        'VariableNames',{'GrowthRate','Frame','Trap'});
end
