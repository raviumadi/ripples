%% field_cr_to_distance_sequences_with_outlier_filter
% -------------------------------------------------------------------------
% PURPOSE
%   Visualise echolocation call timing from field recordings by laying out
%   validated call sequences end-to-end (concatenated across bouts), and
%   converting measured call rates (C_r) into an implied target/cue distance
%   using a responsivity model:
%
%       d = c / [2 (1 + k_r) C_r]
%
%   The script produces a two-panel figure:
%     (1) Call rate (Hz) vs call index (concatenated across sequences)
%     (2) Estimated target distance (m) vs call index
%
%   Points are colour-coded by relative position within each continuous
%   call sequence (start → end), and thin vertical separators mark sequence
%   boundaries.
%
% -------------------------------------------------------------------------
% INPUTS
%   Loads:
%       ../results/DaubAnalysisTable_gcc.mat
%   Expected variable:
%       AnalysisTable  (renamed to T in the script)
%
%   Expected table fields per row (cell arrays):
%       T.Calls{i}            : vector of call IDs / indices for the file
%       T.ValidatedCalls{i}   : subset of Calls that were manually/algorithmically validated
%       T.CallRate{i}         : call rate estimates aligned to ValidatedCalls
%                               (typically length = numValidated - 1)
%       T.CallDurations{i}    : optional; not used in plotting here
%
% -------------------------------------------------------------------------
% SEQUENCE LOGIC (CONTIGUOUS BOUT EXTRACTION)
%   For each table row (recording):
%     1) Map ValidatedCalls into indices within Calls via ismember.
%     2) Identify contiguous runs where successive validated call indices
%        increment by 1 (diff(loc) == 1).
%     3) Treat each contiguous run as one “sequence/bout”.
%     4) For each run, compute call-rate samples and their implied distances.
%
%   Note: Call rate is defined between successive calls, hence for a run of
%   length N calls, the call-rate vector typically has length N-1.
%
% -------------------------------------------------------------------------
% MODEL PARAMETERS
%   c   : speed of sound (m/s), default 343
%   k_r : responsivity factor (dimensionless), default 5
%   fs  : sampling rate (Hz), stored but not used in this plot
%
%   Conversion:
%       callRateToDistance(Cr,kr,c) = c ./ (2 * (1 + kr) .* Cr)
%
% -------------------------------------------------------------------------
% OUTLIER HANDLING
%   Outliers are removed before plotting and summary statistics are computed.
%   Two-layer filtering is applied (see filterOutliers()):
%
%   Rule 1 (robust statistical):
%       - remove outliers in Cr using isoutlier(Cr,'median')
%       - remove outliers in d_est using isoutlier(d_est,'median')
%
%   Rule 2 (hard physical/plot thresholds; configurable):
%       - CR_min < Cr < CR_max
%       - D_min  < d  < D_max
%
%   keepMask = mask1 & mask2 & mask3 & mask4
%
%   The script reports:
%       N      : number of kept points across all sequences
%       N_out  : number removed by filtering
%
% -------------------------------------------------------------------------
% PLOTTING DETAILS
%   Figure: "Sequences laid out sequentially" with two stacked subplots.
%
%   X-axis:
%       A running “call number” index that concatenates all kept points.
%       Each sequence is placed after the previous one with a 1-sample gap.
%
%   Colour:
%       Points are coloured by relative position within each original sequence.
%       (Implemented as a 0→1 ramp; then subset to kept points.)
%
%   Boundaries:
%       Vertical grey lines mark the end of each sequence.
%
%   Colormap:
%       turbo, with a colourbar labelled “Relative position in sequence”.
%
% -------------------------------------------------------------------------
% OUTPUTS
%   1) Figure (not automatically saved unless you enable the saveFigure block)
%   2) Command-window summary statistics across all kept points:
%       Call rate: mean, median, mode, SD, min, max
%       Distance : mean, median, mode, SD, min, max
%
%   3) A stats table is prepared for CSV export:
%       ../fig/field_cr_dist_stats.csv
%     (writetable is currently commented out)
%
% -------------------------------------------------------------------------
% DEPENDENCIES / ASSUMPTIONS
%   - formatLatex(gca) must exist on MATLAB path (plot styling helper).
%   - saveFigure(...) is referenced but commented; if used, it must exist.
%   - The CallRate vectors must be consistent with the contiguous grouping
%     logic (i.e., aligned with ValidatedCalls).
%
% -------------------------------------------------------------------------
% NOTES / INTERPRETATION
%   - The “distance” here is an implied cue distance under the responsivity
%     model; it is not a direct measurement.
%   - Because points are concatenated across many sequences, the x-axis is
%     not time-continuous across sequence boundaries; separators are added
%     to prevent misinterpretation as a single uninterrupted bout.
% -------------------------------------------------------------------------

%% Load data
load('../results/DaubAnalysisTable_gcc.mat'); 
T = AnalysisTable;

% Parameters
c  = 343;        % speed of sound (m/s)
kr = 5;          % responsivity factor
fs = 192e3;      % sampling rate (Hz)

% Outlier thresholds (configurable)
CR_min = 2;      % Hz
CR_max = 200;    % Hz
D_min  = 0;      % m
D_max  = 10;     % m

callRateToDistance = @(Cr,kr,c) c ./ (2 * (1 + kr) .* Cr);

% Prepare figure
figure('Name','Sequences laid out sequentially','Position',[100 100 1000 600]);

ax1 = subplot(2,1,1); hold(ax1,'on'); grid(ax1,'on');
ylabel(ax1,'Call rate (Hz)');
ax1Title = title(ax1,'Call Rates Calculated from Field Recordings');
formatLatex(ax1);

ax2 = subplot(2,1,2); hold(ax2,'on'); grid(ax2,'on');
ylabel(ax2,'Target distance (m)');
ylim(ax2,[0 8])
xlabel(ax2,'Call no. (concatenated across sequences)');
title(ax2,'Estimated Target Distances - Responsivity Model: $d = \frac{c}{2(1+k_r)C_r}, \; k_r = 5$','Interpreter','latex');
formatLatex(ax2);

offset = 0; % keeps track of where to place the next sequence
seqEnds = []; % store sequence end positions for separators
N_total = 0;     % count kept points
N_out_total = 0; % count removed points

% === Store all valid values for stats ===
all_Cr = [];
all_d  = [];

for i = 1:height(T)
    Calls     = T.Calls{i};
    Validated = T.ValidatedCalls{i};
    CallRate  = T.CallRate{i};
    CallDur   = T.CallDurations{i};

    if isempty(Validated) || isempty(CallRate)
        continue;
    end

    [~,loc] = ismember(Validated, Calls);
    loc = loc(loc > 0);
    if isempty(loc), continue; end

    diffs  = diff(loc);
    groups = [0; find(diffs > 1); numel(loc)];

    for g = 1:numel(groups)-1
        idxRange = (groups(g)+1):groups(g+1);
        loc_cont = loc(idxRange);

        if numel(loc_cont) < 2
            continue;
        end

        Cr_valid  = CallRate(idxRange(1:end-1));
        d_est     = callRateToDistance(Cr_valid, kr, c);

        % === Outlier removal ===
        [Cr_filt, d_filt, keepMask] = filterOutliers(Cr_valid, d_est, CR_min, CR_max, D_min, D_max);
        N_total     = N_total + numel(Cr_filt);
        N_out_total = N_out_total + (numel(Cr_valid) - numel(Cr_filt));

        if isempty(Cr_filt), continue; end

        % Collect for stats
        all_Cr = [all_Cr; Cr_filt(:)];
        all_d  = [all_d; d_filt(:)];

        % Gradient colors (aligned with kept points)
        seqColors = linspace(0,1,numel(keepMask))';
        seqColors = seqColors(keepMask);

        x = offset + (1:numel(Cr_filt));  % shift sequence along x

        scatter(ax1, x, Cr_filt, 25, seqColors, 'filled');
        scatter(ax2, x, d_filt,  25, seqColors, 'filled');

        % store sequence end
        seqEnds(end+1) = x(end); %#ok<AGROW>

        offset = offset + numel(Cr_filt) + 1; % add gap between sequences
    end
end

% Add thin vertical lines at sequence boundaries
for s = 1:numel(seqEnds)
    xline(ax1, seqEnds(s)+0.5, '-', 'Color',[0.7 0.7 0.7], 'LineWidth',0.5);
    xline(ax2, seqEnds(s)+0.5, '-', 'Color',[0.7 0.7 0.7], 'LineWidth',0.5);
end

% Update top title with N and N_out
title(ax1, sprintf('Call Rates Calculated from Field Recordings. $N = %d$, $N_{out} = %d$', N_total, N_out_total), ...
    'Interpreter','latex');

% Colormap and colorbar
colormap(turbo)
cb = colorbar(ax2,'eastoutside'); % place outside to the right
cb.Label.String = 'Relative position in sequence (start $\rightarrow$ end)';
cb.Label.Interpreter = 'latex';
cb.TickLabelInterpreter = 'latex';
cb.Box = "off";

% adjust placement further outside
cb.Position(1) = cb.Position(1) + 0.07;
cb.Position(2) = cb.Position(2) + 0.06;
cb.Position(3) = cb.Position(3) * 0.4;   % shrink width to 40% of original
cb.Position(4) = cb.Position(4) * 2;   

%
% fig_path = '../manuscript/fig';
% saveFigure(gcf, fig_path, 'field_cr_dist')

% === Stats ===
cr_mean   = mean(all_Cr);
cr_median = median(all_Cr);
cr_mode   = mode(all_Cr);
cr_std    = std(all_Cr);
cr_min    = min(all_Cr);
cr_max    = max(all_Cr);

d_mean   = mean(all_d);
d_median = median(all_d);
d_mode   = mode(all_d);
d_std    = std(all_d);
d_min    = min(all_d);
d_max    = max(all_d);

fprintf('\n--- Call Rate Stats ---\n');
fprintf('Mean: %.2f Hz | Median: %.2f Hz | Mode: %.2f Hz | SD: %.2f | Min: %.2f | Max: %.2f\n', ...
    cr_mean, cr_median, cr_mode, cr_std, cr_min, cr_max);

fprintf('\n--- Target Distance Stats ---\n');
fprintf('Mean: %.2f m | Median: %.2f m | Mode: %.2f m | SD: %.2f | Min: %.2f | Max: %.2f\n', ...
    d_mean, d_median, d_mode, d_std, d_min, d_max);

% === Write stats to CSV ===
statsTable = table(...
    [cr_mean; d_mean], ...
    [cr_median; d_median], ...
    [cr_mode; d_mode], ...
    [cr_std; d_std], ...
    [cr_min; d_min], ...
    [cr_max; d_max], ...
    'VariableNames', {'Mean','Median','Mode','SD','Min','Max'}, ...
    'RowNames', {'CallRate_Hz','TargetDistance_m'} ...
);

csv_path = '../fig/field_cr_dist_stats.csv';
% writetable(statsTable, csv_path, 'WriteRowNames', true);

% fprintf('\nStats written to %s\n', csv_path);

%% Helper function
function [Cr_filt, d_est_filt, keepMask] = filterOutliers(Cr, d_est, CR_min, CR_max, D_min, D_max)
    % Rule 1: statistical outlier check
    mask1 = ~isoutlier(Cr,'median'); 
    mask2 = ~isoutlier(d_est,'median'); 
    
    % Rule 2: hard thresholds
    mask3 = (Cr > CR_min & Cr < CR_max);
    mask4 = (d_est > D_min & d_est < D_max);
    
    % combine masks
    keepMask = mask1 & mask2 & mask3 & mask4;
    
    % apply
    Cr_filt    = Cr(keepMask);
    d_est_filt = d_est(keepMask);
end