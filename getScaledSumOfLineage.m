function [areas, badSegmentations, detectionTimes,divFrames] = getScaledSumOfLineage(mCell)
% Compute lineage area with division-specific scaling.
%
% The function is similar to getSumOfLineage but a separate scaling factor
% is computed at each cell division. The scaling factor is chosen such that
% the estimated growth rate remains continuous across the division event,
% thereby avoiding artificial jumps in the corrected cell size.
%
% For each division, the scaling factor is computed from the average area
% growth ratios during the two frames preceding the division and the two
% frames following the division:
%
%   scale = mean([A(i-1)/A(i-2), A(i)/A(i-1), ...
%                 D(2)/D(1), D(3)/D(2)]) * A(i) / D(1)
%
% where A denotes the mother-cell area trajectory, D denotes the summed
% daughter-cell area trajectory, and i is the division frame.
%
% Input:
%   mCell - MCell object
%
% Outputs:
%   areas             - Corrected lineage area trajectory.
%   badSegmentations  - Segmentation quality flags corresponding to areas.
%   detectionTimes    - Detection times for each frame.
%   divFrames         - Division-frame indices within the lineage.

if numel(mCell.descendants) == 2
    if mCell.isBadCell == 0
        [a1, bs1, dt1,df1] = getScaledSumOfLineage(mCell.descendants(1));
        [a2, bs2, ~]   = getScaledSumOfLineage(mCell.descendants(2));
        minLength = min(numel(a1), numel(a2));
        if minLength>0
            dareas = a1(1:minLength) + a2(1:minLength);
            % find the scale
            if length(mCell.areas)>=3 && length(dareas)>=3 && all(mCell.badSegmentations(end-2:end)==0) && ...
                    all(bs1(1:3)==0) && all(bs2(1:3)==0)
                delta_array = [mCell.areas(end-1:end)./mCell.areas(end-2:end-1) ...
                             dareas(2:3)./dareas(1:2)];
                scale = mean(delta_array) * mCell.areas(end) / dareas(1);
                dareas = scale * dareas;
            end
            dbs = bs1(1:minLength) + bs2(1:minLength);
            ddt = dt1(1:minLength);            
            areas = [mCell.areas, dareas];
            badSegmentations = [mCell.badSegmentations, dbs];
            detectionTimes = [mCell.detectionTimes; ddt];
            divFrames = [mCell.lifeTime mCell.lifeTime+df1];
            return;
        end
    else
        areas = [];
        badSegmentations = [];
        detectionTimes = [];
        divFrames = [];
        return
    end
end

if mCell.isBadCell == 0 %|| mCell.isBadCell == 4 %|| mCell.isBadCell == 3
    areas = mCell.areas;
    badSegmentations = mCell.badSegmentations;
    detectionTimes = mCell.detectionTimes;
    divFrames = mCell.lifeTime;
elseif mCell.isBadCell == 2
    stop = find(mCell.badSegmentations==2,1);
    areas = mCell.areas(1:stop);
    badSegmentations = mCell.badSegmentations(1:stop);
    detectionTimes = mCell.detectionTimes(1:stop);
    divFrames = stop;
else
    areas = [];
    badSegmentations = [];
    detectionTimes = [];
    divFrames = [];
end