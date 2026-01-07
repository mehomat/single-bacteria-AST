function Im = getKymograph(eln, sheet, positionRes, trapRes, positionSus, trapSus, endFrame, kymographfilename, varargin)

% getKymograph
% Kymograph montage
%
% Input arguments:
% - eln: string/char, microfluidic ELN
% - sheet: string/char, sheet name in excel
% - positionRes: scalar, position for resistant strain
% - trapRes: scalar, trap for resistant strain
% - positionSus: scalar, position for susceptible strain
% - trapSus: scalar, trap for susceptible strain
% - endFrame: scalar, end frame of the experiment (must be divisible by 5)
% - kymographfilename: string/char, output tiff filename to write
%
% Name-value pairs:
% - 'CodeDir': string/char, default 'code'
% - 'OutputDir': string/char, required by the wrapper
%
% Output:
% T: Kymograph montage

%% Parse inputs

p = inputParser;
p.addParameter('CodeDir', 'code');
p.addParameter('OutputDir', '');
p.parse(varargin{:});

outputDir = p.Results.OutputDir;

%% Load experiment data from OutputDir

expInfoObj = loadExpInfo(outputDir);

%% Susceptible strain kymograph
p = positionSus;
trap = trapSus;
imSus = displayCellsInMMTrapPhase(expInfoObj, p, trap,'Range', 1:5:endFrame, 'ShowOutlines', 0);

%% Resistant strain kymograph
p = positionRes;
trap = trapRes;
imRes = displayCellsInMMTrapPhase(expInfoObj, p, trap, 'Range', 1:5:endFrame, 'ShowOutlines', 0);

%% Rotate images that are not correct

incorrect = ["EXP-25-CB7064"];

if ismember(eln, incorrect)

    if eln == incorrect(1)

        imSus = flipud(imSus);
        imRes = flipud(imRes);
    
    end

end

%% Crop
[rows1, cols1] = size(imRes);
[rows2, cols2] = size(imSus);
targetCols = min(cols1, cols2);
imRes = imRes(:, 1:targetCols);
imSus = imSus(:, 1:targetCols);

%% Choose rows to extract
rows_res = max(1,700):min(rows1,1470);
rows_sus = max(1,700):min(rows2,1470);

% Contrast adjust
imResCrop = imadjust(imRes(rows_res, :));
imSusCrop = imadjust(imSus(rows_sus, :));

%% Make montage

if sheet == "Cell Elongation"
    Im = [imSusCrop; imResCrop];
else
    Im = [imResCrop; imSusCrop];
end

%% Save montage
imwrite(Im, kymographfilename);

%% Display
figure;
imshow(Im);
end
