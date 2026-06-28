function T = computeDivisionSizes(outputDir, cutoffFrame, posIdxRange, varargin)

% computeDivisionSizes
% For N = 2:4, finds all valid lineages, then stores [Gen1, GenLast, Gap]
%
% Input:
% - outputDir: string/char, where the expInfoObj is loaded from
% - cutoffFrame: scalar, exclude cells with lastFrame > cutoffFrame
% - posIdxRange: numeric indices into expInfoObj.getPositionList()
%
% Name-value pairs:
% - lengthCutOff: cell tracks
% - parentLengthCutOff: parent length
% - daughterLengthCutOff: daughter length
%
% Output
% - T: table with columns: Gen 1, Gen 2, Gap, Position
% (Gen 2 means "last generation". Meaning Gen2/Gen3/Gen4 depending on Gap)

    % Parse optional cutoffs
    lengthCutOff = 10;
    parentLengthCutOff = 3;
    daughterLengthCutOff = 3;

    if numel(varargin) >= 1, lengthCutOff = varargin{1}; end
    if numel(varargin) >= 2, parentLengthCutOff = varargin{2}; end
    if numel(varargin) >= 3, daughterLengthCutOff = varargin{3}; end

    expInfoObj = loadExpInfo(outputDir);
    posList = expInfoObj.getPositionList();

    if ~isempty(posIdxRange)
        posList = posList(posIdxRange);
    end

    allGen1 = [];
    allGen2 = [];
    allGap = [];
    allPos = strings(0,1);

    for pIdx = 1:numel(posList)
        pos = posList{pIdx};

        trackedFile = fullfile(outputDir, pos, "trackedCells.mat");
        if ~isfile(trackedFile)
            continue;
        end

        mCells = Cell.MCell.loadMCells(trackedFile);
        if isempty(mCells)
            continue;
        end

        lastFrames = [mCells.lastFrame];
        [divAreas, genIndex, divFrames] = getDivisionAreasForDivSizes(mCells, lengthCutOff, daughterLengthCutOff);

        divAreas2 = nan(size(lastFrames));
        divFrames2 = nan(size(lastFrames));

        divAreas2(genIndex) = divAreas;
        divFrames2(genIndex) = divFrames;

        % Exclude division events post cutoff
        divAreas2(divFrames2 > cutoffFrame) = nan;

        for N = 2:4
            lineages = selectCellLineages(mCells, N, lengthCutOff, parentLengthCutOff, daughterLengthCutOff);

            % Added for verification
            fprintf('Position %s | N=%d | Raw lineages: %d\n', ...
                pos, N, size(lineages,1));

            if isempty(lineages)
                continue;
            end

            linDivAreas = nan(size(lineages)); % lineageCount x N
            for i = 1:size(lineages, 1)
                c = mCells(lineages(i,1));
                linDivAreas(i,1) = divAreas2(c.id);

                for j = 2:N
                    c = c.descendants(lineages(i,j));
                    linDivAreas(i,j) = divAreas2(c.id);
                end
            end

            % Keep only full valid lineages
            good = ~any(isnan(linDivAreas), 2);
            linDivAreas = linDivAreas(good, :);

            % Added for verification
            fprintf('Position %s | N = %d | Valid after area filter: %d\n', ...
                pos, N, size(linDivAreas,1));
            if isempty(linDivAreas)
                fprintf('Position %s | N = %d | No surviving full lineages\n', pos, N);
                continue;
            end 

            % Keep first and last only, add gap
            g1 = linDivAreas(:, 1);
            g2 = linDivAreas(:, end);
            gap = (N-1) * ones(size(g1));

            allGen1 = [allGen1; g1];
            allGen2 = [allGen2; g2];
            allGap = [allGap; gap];
            allPos = [allPos; repmat(string(pos), numel(g1), 1)];
        end
    end

    T = table(allGen1, allGen2, allGap, allPos, ...
        'VariableNames', {'Gen 1','Gen 2','Gap','Position'});
end
