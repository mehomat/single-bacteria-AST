function T = computeGenerationTimes(outputDir, endFrame, posIdxRange, lengthCutOff, parentLengthCutOff, daughterLengthCutOff)

% Per cell generation times 
%
% Loops over positions, keeps only cells whose track
% ends on or before endFrames and measures the generation time of each cell
% with getGenerationTimes
%
% Input:
% - outputDir: path to the experiment's output folder 
% - endFrame: scalar frame number, cells whose last frame is after 
% endFrame are excluded
% - posIdxRange: vector of position indices to loop over
% - lengthCutOff: minimum cell track length to include (default 5)
% - parentLengthCutOff: minimum parent track length (default 3)
% - daughterLengthCutOff: minimum daughter track length (default 3)
%
% Output:
% - T: table with one row per cell        

    if nargin < 4 || isempty(lengthCutOff); lengthCutOff = 5; end
    if nargin < 5 || isempty(parentLengthCutOff); parentLengthCutOff = 3; end
    if nargin < 6 || isempty(daughterLengthCutOff); daughterLengthCutOff = 3; end

    expInfoObj = loadExpInfo(outputDir);

    pos_all = [];
    idx_all = [];
    gt_all = [];

    for pos = posIdxRange
        mCells = expInfoObj.getMCells(pos);
        if isempty(mCells); continue; end

        % Exclude any cell whose last frame is after endFrame
        keep = arrayfun(@(c) ~isempty(c.lastFrame) && c.lastFrame <= endFrame, mCells);
        mCells = mCells(keep);
        if isempty(mCells); continue; end

        [gt, cellIdx] = getGenerationTimes(mCells, lengthCutOff, parentLengthCutOff, daughterLengthCutOff);

        gt_all = [gt_all, gt];                       
        idx_all = [idx_all, cellIdx];                   
        pos_all = [pos_all, repmat(pos, 1, numel(gt))]; 

    end

    T = table(pos_all(:), idx_all(:), gt_all(:), ...
        'VariableNames', {'Position', 'CellIndex', 'GenerationTime_min'});
end