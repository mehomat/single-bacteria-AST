function [areas, badSegmentations, detectionTimes,divFrames] = getSumOfLineage(mCell)

if numel(mCell.descendants) == 2
    if mCell.isBadCell == 0
        [a1, bs1, dt1,df1] = getSumOfLineage(mCell.descendants(1));
        [a2, bs2, ~]   = getSumOfLineage(mCell.descendants(2));
        minLength = min(numel(a1), numel(a2));
        if minLength>0
            dareas = a1(1:minLength) + a2(1:minLength);
            dbs = bs1(1:minLength) + bs2(1:minLength);
            ddt = dt1(1:minLength);            
            areas = [mCell.areas, dareas];
            %areas = [mCell.lengths, dareas];
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
    %areas = mCell.lengths(1:stop);
    badSegmentations = mCell.badSegmentations(1:stop);
    detectionTimes = mCell.detectionTimes(1:stop);
    divFrames = stop;
else
    areas = [];
    badSegmentations = [];
    detectionTimes = [];
    divFrames = [];
end