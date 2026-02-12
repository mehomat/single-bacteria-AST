function [areas, cellIndices, divFrames] = getDivisionAreasForDivSizes(mCells,varargin)

% This function extracts information about cell division areas. The cell
% should have two daughter cells, all the cell tracks
% should be sufficiently long. Then the division area is the area of the 
% cell in the last detection frame.
%
% Input:
% - mCells: array of MCell objects
% - lengthCutOff: cutoff threshold of cell track length, optional, 5 by default
% - daughterLengthCutoff: cutoff threshold of daughter cell track length,
%   optional, 3 by default
%
% getDivisionAreas(...,'PixelSize',pixelSize) scales the output values from
% px^2 to um^2 according to pixelSize (in um). Optional, no scaling by
% default (pixelSize=1).
%
% Output:
% - areas: array with cell division areas
% - cellIndices: array with indices of the cells where division area was measured
% - divFrames: array with last frames

ip = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
ip.addOptional('lengthCutOff', 5, validScalarPosNum);
ip.addOptional('daughterLengthCutOff',3,validScalarPosNum);
ip.addParameter('PixelSize',1,validScalarPosNum)
ip.parse(varargin{:});
lengthCutOff = ip.Results.lengthCutOff;
daughterLengthCutOff = ip.Results.daughterLengthCutOff;
pixelSize = ip.Results.PixelSize;

areas = nan(1, length(mCells));
divFrames = nan(1, length(mCells)); % Initialize storage of division frames

for i = 1:length(mCells)
    c = mCells(i);
     if length(c.descendants) == 2 && all([c.descendants.lifeTime] >= daughterLengthCutOff)
        if c.isBadCell == 0 && c.lifeTime >= lengthCutOff
            if c.badSegmentations(end) == 0
                areas(i) = c.areas(end);
                divFrames(i) = c.lastFrame; 
            end     
        end
     end
end

cellIndices = find(~isnan(areas));
areas = areas(cellIndices);
divFrames = divFrames(cellIndices);

if pixelSize~=1
    areas = areas * pixelSize^2;
end