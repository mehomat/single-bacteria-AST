function T = getTrapLengths(expInfoObj,switchFrame,strain,posRange,yThresh)
% T = getTrapLengths(expInfoObj,switchFrame,strain,posRange,yThresh)
% Computes the length of bacterial cell occupancy within each mother-
% machine trap over time.
%
% Input:
%   expInfoObj  - expInfo object containing image analysis information
%   switchFrame - Frame index corresponding to the first media switch.
%   strain      - (Optional) String specifying the strain name. Default: ""
%   posRange    - (Optional) Array specifying which positions to include. 
%                 Default: include all positions.
%   yThresh     - (Optional) A y-axis threshold used to discard trap rows 
%                 beyond a specified region. Default: [] (no thresholding).
%
% Output:
%   T - Table with the following columns:
%         'Length'       : Measured trap length (in pixels).
%         'Frame'        : Frame index.
%         'Trap'         : Unique trap identifier, defined as
%                          (position - 1) * nTraps + trapIndex.
%         'CellsAtStart' : Number of cells at the start frame
%         'CellsAtSwitch': Number of cells at the switch frame
%         'Strain'       : Strain name.

if nargin<3, strain = "";end
if nargin<4, posRange = [];end
if nargin<5, yThresh = [];end

posList = expInfoObj.getPositionList();
if isempty(posRange), posRange=1:length(posList);end

outputDir = expInfoObj.getPathName('output');
param = expInfoObj.getParameters();
imSize = param.roiInPhase([4 3])+1;
nTraps = param.nGrowthChannels-length(param.emptyChannel);
nFrames = min(cellfun("length",expInfoObj.imRange));
nPositions = length(posRange);

allTrapLength = cell(1,nPositions);
cellCountsAtStart = zeros(nTraps,nPositions);
cellCountsAtSwitch = zeros(nTraps,nPositions);

parfor pi=1:length(posRange)
    pos = posRange(pi);
    fprintf('Position %d/%d\n',pos,nPositions);
    posName = posList{pos};
    traps = expInfoObj.getChannelLocations(posName);
    if ~isempty(param.emptyChannel)
        trapIndices = 1:param.nGrowthChannels;
        trapIndices(param.emptyChannel(2:end)) = [];
        emptyChIdx = find(trapIndices == param.emptyChannel(1));
        traps = traps([1:emptyChIdx-1 emptyChIdx+1:end],:,:);
    end
    trapMask = zeros(imSize);
    for t=1:nTraps
        trapMask(traps(t,2):traps(t,2)+traps(t,4),traps(t,1):traps(t,1)+traps(t,3))=t;
    end
    if ~isempty(yThresh)
        if yThresh>imSize(1)/2
            trapMask(yThresh:end,:)=0;
        else
            trapMask(1:yThresh,:)=0;
        end
    end
    segDir = fullfile(outputDir,posName,'SegmentedPhase');
    segList = dir(fullfile(segDir,'*.tif*'));
    trapLengths = zeros(nFrames,nTraps);
    for i=1:nFrames
        im = imread(fullfile(segDir,segList(i).name));
        imTrapMask = trapMask.*(im>0);
        props = regionprops(imTrapMask,'BoundingBox');
        for t=1:nTraps
            if t<=length(props)
                bb = props(t).BoundingBox;
                if ~isempty(bb)
                    trapLengths(i,t)=bb(4);
                end
            end
        end
    end
    allTrapLength{pi}=trapLengths;

    % compute cell counts
    mCells = expInfoObj.getMCells(pi);
    startFrames = [mCells.birthFrame];
    lastFrames = [mCells.lastFrame];
    growthChannels = [mCells.growthChannel];
    isAtStartFrame = startFrames==1;
    isAtSwitchFrame = startFrames<=switchFrame & lastFrames>=switchFrame;
    for t=1:nTraps
        cellCountsAtStart(t,pi)=sum(isAtStartFrame & growthChannels==t);
        cellCountsAtSwitch(t,pi)=sum(isAtSwitchFrame & growthChannels==t);
    end
end

% Convert cell arrays to a matrix
M = zeros(nTraps*nPositions*nFrames,5);
frames = (1:nFrames)';
ff = ones(nFrames,1);
for pi=1:length(posRange)
    pos = posRange(pi);
    for t=1:nTraps
        trapID = (pos-1)*nTraps+t;
        trapID2 = (pi-1)*nTraps+t;
        ind = nFrames*(trapID2-1)+1:nFrames*trapID2;
        M(ind,:) = [allTrapLength{pi}(:,t) frames trapID*ff, ...
                    cellCountsAtStart(t,pi)*ff cellCountsAtSwitch(t,pi)*ff];
    end
end

% Convert matrix to a table
T = array2table(M,"VariableNames",["Length","Frame","Trap","CellsAtStart","CellsAtSwitch"]);
T = [T array2table(repmat(strain,size(T,1),1),"VariableNames","Strain")];