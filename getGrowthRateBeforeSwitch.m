function [gr, grStd] = getGrowthRateBeforeSwitch(mCell, switchFrame,makePlot)
% Fit single-term exponential function 
%   y = a * exp(b*t)
% to find cell area growth rate before media switch.
% Parameter b determines the growth rate.
% 
% The cell track must have "OK" segmentation labels in 
% at least 75% of frames.
%
% Input:
%   mCell - MCell object with cell track info
%   switchFrame - the frame when media is switched to, e.g., antibiotics.
%   makePlot - logical, if 1 plot fitted curve and input data. 
%              Optional, 0 by default.
%
% Output:
%   gr - fitted growth rate
%   grStd - standard deviation from 66.(6)% confidence interval

if nargin<3
    makePlot = 0;
end

selFrames = 1:min(mCell.lifeTime,switchFrame-mCell.birthFrame+1);
badSegs = mCell.badSegmentations(selFrames);
if switchFrame<mCell.birthFrame+3 || sum(~badSegs) < max(3,0.75*length(badSegs))
    gr = nan;
    grStd = nan;
    return
end
areas = mCell.areas(selFrames);
detectionTimes = mCell.detectionTimes(selFrames);
relDetectionTimes =  detectionTimes - detectionTimes(1);
relDetectionTimes = relDetectionTimes / (60 * 1000); % convert to minutes

y = areas(badSegs == 0);
t = relDetectionTimes(badSegs==0);

fitObj = fit(t,y', 'exp1');
gr = fitObj.b;

if nargout == 2
    intervals = confint(fitObj, 0.66666666);
    devFromMean = intervals(:,2) - fitObj.b;
    grStd = devFromMean(2);
end

if makePlot
    figure,
    plot(fitObj,t,y)
end