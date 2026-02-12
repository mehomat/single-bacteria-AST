function lineages = selectCellLineages(mCells, N, varargin)

% Searches for cell lineages that links N cells
%
% Each cell in the lineage must satisfy the following conditions:
% 1) to be an "OK" cell.
% 2) Parent rule:
% - If the cell has a parent and parentLengthCutOff > 0, the parent must be "OK".
% - If the cell has no parent (founder cell at experiment start), it is still allowed.
% 3) If daughterLengthCutOff > 0, the cell has two "OK" daughter cells.
% 4) If sisterLengthCutOff > 0, the cell has an "OK" sister cell.
%
% A cell track is considered "OK" if
% - it has labelled isBadCell = 0;
% - its length >= given length cutoff
% - it has labels badSegmentations=0 in >= 90 % frames.
%
% Input:
% - mCells: array of MCell objects
% - N: number of cell generations in a lineage
%
% Name-value pairs:
% - lengthCutOff: cutoff threshold of cell track length, 5 by default
% - parentLengthCutoff: cutoff threshold of parent cell track length, 3 by default
% - daughterLengthCutoff: cutoff threshold of daughter cell track length, 3 by default
% - sisterLengthCutoff: cutoff threshold of sister cell track length, 0 by default
% - useDaughters: if 1, include to the lineage daughter cells of the Nth
% cell in the lineage. Could be handy for computing single
% cell initiation area. the default value is 0.
% 
% Output:
% - lineages: matrix of size M x K, each row represents a cell lineage:
% cell id, index of its daughter cell, index of its daughter's
% daughter cell etc. The indices are either 1 or 2.
% K = N if useDaughters=0 and K = N+1 otherwise.

    % Parse parameters
    ip = inputParser;

    validScalarNonNegNum = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
    ip.addOptional('lengthCutOff',5,validScalarNonNegNum);
    ip.addOptional('parentLengthCutOff',3,validScalarNonNegNum);
    ip.addOptional('daughterLengthCutOff',3,validScalarNonNegNum);
    ip.addOptional('sisterLengthCutOff',0,validScalarNonNegNum);
    ip.addOptional('useDaughters',0,validScalarNonNegNum);
    ip.parse(varargin{:});

    L = ip.Results.lengthCutOff;
    pL = ip.Results.parentLengthCutOff;
    dL = ip.Results.daughterLengthCutOff;
    sL = ip.Results.sisterLengthCutOff;
    useDaughters = ip.Results.useDaughters;

    lineages = [];
    for i = 1:length(mCells)
        c = mCells(i);
        lineage = selectLineage(c, N, L, pL, dL, sL, 0);
        if ~isempty(lineage)
            lineages = [lineages; lineage];
        end
    end

    if useDaughters
        tmpLineages = ones(2*size(lineages, 1), N+1);
        tmpLineages(1:2:end, 1:N) = lineages;
        tmpLineages(2:2:end, 1:N) = lineages;
        tmpLineages(2:2:end, N+1) = 2;
        lineages = tmpLineages;
    end
end


% Utility functions
function lin = selectLineage(c, N, L, pL, dL, sL, branch)
    lin = [];
    if isOK(c,L,1)
        if isempty(c.parent) || isOK(c.parent, pL)
            if ~isempty(c.parent)
                if c.parent.descendants(1).id == c.id
                    sc = c.parent.descendants(2);
                else
                    sc  = c.parent.descendants(1);
                end
            else
                sc = [];
            end
            if isOK(sc, sL)
                dcs = c.descendants;
                if N > 1 && length(dcs) == 2
                    if isOK(dcs,dL)
                        for i = 1:2
                            d_lin = selectLineage(dcs(i), N-1, L, pL, dL, sL, i);
                            if ~isempty(d_lin)
                                if branch > 0
                                    col = repmat(branch,size(d_lin,1),1);
                                else
                                    col = repmat(c.id,size(d_lin,1),1);
                                end
                                lin = [lin; col d_lin];
                            end
                        end
                    end
                elseif N==1 && isOK(dcs, dL)
                    if branch > 0
                        lin = branch;
                    else
                        lin = c.id;
                    end
                end
            end
        end
    end
end

function status = isOK(mCells, lengthCutOff,checkBadCell)
    if nargin < 3
        checkBadCell = 0;
    end
    if ~lengthCutOff
        status = true;
    elseif ~isempty(mCells)
        status = true;
        for i=1:length(mCells)
            c = mCells(i);
            if  (checkBadCell && c.isBadCell>0) || c.lifeTime<lengthCutOff || ...
                    sum(c.badSegmentations==0)<0.9*c.lifeTime
                status = false;
                return
            end
        end
    else
        status = false;
    end
end