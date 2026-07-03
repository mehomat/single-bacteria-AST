function [grs,frs] = getSingleLineageGrowthRate(mCell,window,dt,fittype,previousLastFrame)
% Calculates growth rates for a cell lineage and, recursively, for one of
% its daughter lineages. Growth rates are estimated by sliding-window 
% exponential fits to the summed segmentation areas of the lineage.
% When cell divides, the sum of daughter cell areas is scaled to match the
% average growth rate in the 2 frames before the division and 2 frames
% after the division.
% The growth-rate trace continues through the daughter with the longer 
% lineage. A minimum of 5 frames is required for a fit.
%
% INPUTS
%   mCell             - MCell object representing a tracked cell.
%   window            - 1×2 vector [wPre, wPost] defining the half-widths
%                       of the sliding window in frames.
%   dt                - Scalar, time interval between frames in minutes.
%   fittype           - Fit type: 'exp1' (exponential) or 'poly1' (linear)
%   previousLastFrame - (Optional) Frame index of the last frame used
%                       by the preceding lineage segment to avoid overlap
%                       in growth rate computation. Default: [].
%
% OUTPUTS
%   grs  - array of instantaneous growth rates.
%   frs  - array of frame indices corresponding to each entry in grs.
%
% See also: getScaledSumOfLineage, movgrowthrate2.

if nargin<5, previousLastFrame = [];end

minLength = 5;  % minimum number of frames required to attempt a fit

% Extract area time series
[areas, badSegmentations] = getScaledSumOfLineage(mCell);
frs = mCell.birthFrame : (mCell.birthFrame + length(areas) - 1);

% Remove frames that fall within window(1) frames of previousLastFrame to
% avoid fitting the same data twice.
if ~isempty(previousLastFrame)
    firstFrameIndex = find(frs > previousLastFrame - window(1), 1);
    if ~isempty(firstFrameIndex)
        areas = areas(firstFrameIndex:end);
        badSegmentations = badSegmentations(firstFrameIndex:end);
        frs = frs(firstFrameIndex:end); 
    else
        frs = [];
    end    
end

% Fit growth rates
if numel(frs) > minLength
    t = frs * dt;
    grs = movgrowthrate2(t, areas, badSegmentations, window, fittype);

    % Retain only frames strictly after previousLastFrame
    if ~isempty(previousLastFrame)
        firstFrameIndex = find(frs > previousLastFrame, 1);
        grs = grs(firstFrameIndex:end);
        frs = frs(firstFrameIndex:end);
    end 

    % Clean up
    firstFiniteIndex = find(isfinite(grs),1);
    lastFiniteIndex = find(isfinite(grs),1,'last');
    if isempty(firstFiniteIndex)
        indFinite = [];
    else
        indFinite = firstFiniteIndex:lastFiniteIndex;
    end
    grs = grs(indFinite);
    frs = frs(indFinite);
else
    frs = [];
    grs = [];
    return;
end

% Ensure column vectors
frs = frs(:);
grs = grs(:);

% If the cell divided, continue the growth-rate trace through the daughter
% whose lineage is longer
if ~isempty(frs) && numel(mCell.descendants) == 2
    a1 = getSumOfLineage(mCell.descendants(1));
    a2 = getSumOfLineage(mCell.descendants(2));
    if numel(a1)>numel(a2)
        daughterCell = mCell.descendants(1);
    else
        daughterCell = mCell.descendants(2);
    end
    [dgrs,dfrs] = getSingleLineageGrowthRate(daughterCell, window, dt, ....
                                           fittype, frs(end));
    frs = [frs; dfrs];
    grs = [grs; dgrs];
end