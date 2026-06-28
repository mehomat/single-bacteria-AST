function T = makeTrajectoryLossTable(csvDir, varargin)

% All per ELN trajectory loss summaries into one table
%
% Reads every file matching the pattern in csvDir (one row per ELN, as
% written by getPreSwitchTrajectoryLoss), concatenates them, sorts by ELN,
% optionally appends a pooled TOTAL row (counts summed, percentages
% recomputed) and writes a combined CSV.
%
% Input:
% - csvDir: folder containing the per-ELN summary CSVs
%
% Name-value pairs:
% - Pattern: filename glob, default "trajectory_loss_summary_*.csv"
%            (the ALL file is excluded automatically)
% - OutFile: path for the combined CSV, default
%            fullfile(csvDir,"trajectory_loss_summary_ALL.csv")
% - AddTotal:logical, append a pooled TOTAL row, default true
%
% Output:
% - T: combined table (also written to OutFile)

p = inputParser;
p.addParameter('Pattern', "trajectory_loss_summary_*.csv");
p.addParameter('OutFile', "");
p.addParameter('AddTotal', true, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
pattern  = string(p.Results.Pattern);
outFile  = string(p.Results.OutFile);
addTotal = logical(p.Results.AddTotal);

if outFile == ""
    outFile = fullfile(csvDir, "trajectory_loss_summary_ALL.csv");
end

files = dir(fullfile(csvDir, pattern));
if isempty(files)
    error('makeTrajectoryLossTable: no files matching %s in %s', pattern, csvDir);
end

%% Read and concatenate (skip the combined ALL file if it is matched)
T = table();
for i = 1:numel(files)
    if contains(files(i).name, "_ALL.csv")
        continue;
    end
    Ti = readtable(fullfile(files(i).folder, files(i).name), 'TextType', 'string');
    T = [T; Ti];
end

if height(T) == 0
    error('makeTrajectoryLossTable: no per-ELN summaries read.');
end

% Sort by ELN 
if ismember('ELN', T.Properties.VariableNames)
    T = sortrows(T, 'ELN');
end

%% Optional pooled total row (sum counts, recompute percentages)
if addTotal
    tot = table();
    if ismember('ELN', T.Properties.VariableNames)
        tot.ELN = "TOTAL";
    end
    tot.SwitchFrame = NaN; % not meaningful pooled
    tot.PreSwitchWindowMin = T.PreSwitchWindowMin(1); % constant across ELNs
    tot.nTrapsTotal = sum(T.nTrapsTotal);
    tot.nNonEmptyTraps = sum(T.nNonEmptyTraps);
    tot.nValidTraps = sum(T.nValidTraps);
    tot.pctValidOfNonEmpty = 100 * tot.nValidTraps / tot.nNonEmptyTraps;
    tot.nValidLineages = sum(T.nValidLineages);
    tot.nLoss = sum(T.nLoss);
    tot.LossRate_pct = 100 * tot.nLoss / tot.nValidLineages;
    tot.nLostBeforeSwitch = sum(T.nLostBeforeSwitch);
    tot.pctLostBeforeSwitch = 100 * tot.nLostBeforeSwitch / tot.nValidTraps;
    tot.nLostBeforeABEffect = sum(T.nLostBeforeABEffect);
    tot.pctLostBeforeABEffect = 100 * tot.nLostBeforeABEffect / tot.nValidTraps;
    tot.nNoData = sum(T.nNoData);
    tot.nFilled = sum(T.nFilled);
    tot.nFullLength = sum(T.nFullLength);
    tot.nFailsGenTime = sum(T.nFailsGenTime);
    tot.pctFailsGenTime = 100 * tot.nFailsGenTime / tot.nValidTraps;
    tot.nGenTimeFailButComplete = sum(T.nGenTimeFailButComplete);
    tot.nFailsBaseline = sum(T.nFailsBaseline);
    tot.nInFinalAnalysis = sum(T.nInFinalAnalysis);

    % Align column order to T and append
    tot = tot(:, T.Properties.VariableNames);
    T = [T; tot];
end

%% Write
writetable(T, outFile);
fprintf('makeTrajectoryLossTable: wrote %d rows -> %s\n', height(T), outFile);

end