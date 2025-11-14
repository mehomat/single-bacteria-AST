function T = getSingleCellLineageAreas(mCell,strain)
% T = getSingleCellLineageAreas(mCell,strain)
% Computes the cell area for a single lineage over time. 
% When a division occurs, the function selects the daughter cell that is 
% closest to the mother-machine trap constriction.
%
% Input:
%   mCell  - MCell object containing cell lineage information.
%   strain - (Optional) String specifying the strain name. Default: "".
%
% Output:
%   T - A table with the following columns:
%         'Frame'  : Frame index.
%         'Area'   : Measured trap length (in pixels).
%         'Strain' : Strain name.

if nargin<2
    strain="";
end

% detect cell growth orientation along y-axis.
signy=sign(mean(diff(mCell.centroids(mCell.badSegmentations==0,2))));  

areas=[];
frames = [];
while ~isempty(mCell) && (mCell.isBadCell==0 || mCell.isBadCell==2)% || mCell.isBadCell==3)
    if mCell.isBadCell==2
        ind = 1:find(mCell.badSegmentations>0,1)-1;
    else
        ind = mCell.badSegmentations==0;
    end
    cframes = mCell.birthFrame:mCell.lastFrame;
    frames=[frames cframes(ind)];
    areas = [areas mCell.areas(ind)];
    if isempty(mCell.descendants)
        mCell=[];
    else
        dc1= mCell.descendants(1);
        idx1=find(dc1.badSegmentations==0,1);
        if ~isempty(idx1)
            dc1_y=dc1.centroids(idx1,2);
        else
            dc1_y=nan;
        end
        dc2= mCell.descendants(2);
        idx2=find(dc2.badSegmentations==0,1);
        if ~isempty(idx2)
            dc2_y=dc2.centroids(idx2,2);
        else
            dc2_y=nan;
        end
        if all(isnan([dc1_y,dc2_y]))
            break
        else
            if signy*dc1_y < signy*dc2_y
                mCell=dc1;
            else
                mCell=dc2;
            end
        end
    end
end

T=table(frames',areas',repmat(strain,size(frames')),'VariableNames',["Frame","Area","Strain"]);
