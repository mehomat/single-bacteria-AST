function Im = getKymographForFig1(pos, trap, endFrame, upperLimit, lowerLimit, kymographfilename, outputDir)

% getKymographForFig1
% Create kymograph montage for one position/trap:
% [phase kymograph, segmented mask kymograph]
%
% Inputs:
% - pos: scalar position (the actual position number to read from expInfo)
% - trap: scalar trap number
% - endFrame: scalar end frame of experiment
% - upperLimit: scalar row index
% - lowerLimit: scalar row index
% - kymographfilename: string/char
% - outputDir: string/char, where the expInfoObj is loaded from
%
% Output:
% - Kymograph montage including phase contrast image and segmented mask

% Load expInfoObj
expInfoObj = loadExpInfo(outputDir);

imPhase = displayCellsInMMTrapPhase(expInfoObj, pos, trap, ...
    'Range', 1:5:endFrame, 'ShowOutlines', 0);

imSeg = displayCellsInMMTrapPhase(expInfoObj, pos, trap, ...
    'Range', 1:5:endFrame, 'ShowSegmentation', 1, 'ShowOutlines', 0);


% Sizes 
rows1 = size(imPhase, 1); cols1 = size(imPhase, 2);
rows2 = size(imSeg, 1); cols2 = size(imSeg, 2);

% Match columns 
targetCols = min(cols1, cols2);
imPhase = imPhase(:, 1:targetCols);
imSeg = imSeg(:, 1:targetCols);

% Crop 
upper = max(1, upperLimit);
lower1 = min(rows1, lowerLimit);
lower2 = min(rows2, lowerLimit);

imPhaseCrop = imadjust(imPhase(upper:lower1, :));
imSegCrop = imSeg(upper:lower2, :);

Im = [imPhaseCrop; imSegCrop];

imwrite(Im, kymographfilename);

end

