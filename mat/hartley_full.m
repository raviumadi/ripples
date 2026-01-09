%% hartley_fig9_responsivity_fit_from_CR_then_invert_IPI_UNION
% -------------------------------------------------------------------------
% PURPOSE
%   Reanalyse Hartley et al. (1989) call-timing data using a simple
%   responsivity model, while UNION-augmenting the datasets so that:
%     (i)  all call-rate (CR) points are shown, even where no IPI exists, and
%     (ii) all inter-pulse-interval (IPI) points are shown, even where no CR exists.
%   Missing “counterpart” values (IPI at CR distances; CR at IPI distances)
%   are filled by interpolation versus distance, with a conservative gap
%   limit to avoid uncontrolled extrapolation.
%
%   Output figures:
%     Figure 1: CR vs reported distance (measured CR + filled CR-at-IPI-dist)
%     Figure 2: IPI-implied cue distance vs reported distance
%               (measured IPI inversion + filled IPI-at-CR-dist inversion)
%
% -------------------------------------------------------------------------
% BACKGROUND / MODEL
%   Let c be speed of sound. Under a responsivity model (gain-like parameter k_r),
%   the “effective” cue distance inferred from IPI is:
%
%       d_est = (c * IPI) / (2 * (1 + k_r))
%
%   and the predicted call rate as a function of distance is:
%
%       CR_pred(d) = c / (2 * (1 + k_r) * d)
%
%   Here, k_r is anchored from the high-rate / near-target regime of the CR data.
%
% -------------------------------------------------------------------------
% DATA INPUTS
%   Requires digitised CSV files:
%     - CR vs distance:
%         data/hartley_1989/Bat1_CR_Distance.csv
%         data/hartley_1989/Bat2_CR_Distance.csv
%     - IPI vs distance:
%         data/hartley_1989/Bat_1_IPI_Distance.csv
%         data/hartley_1989/Bat2_IPI_Distance.csv
%
%   Each CSV must have at least two columns:
%     col1 = distance (m)
%     col2 = CR (Hz)   OR  IPI (ms)
%
% -------------------------------------------------------------------------
% UNION AUGMENTATION / “FILL” STRATEGY
%   Two complementary fills are performed using interpolation vs distance:
%
%   A) For every measured CR point at distance d_CR:
%        - Fill IPI(d_CR) by interpolating IPI vs distance (Bat-specific)
%        - Convert to d_est_fill using the anchored k_r
%
%   B) For every measured IPI point at distance d_IPI:
%        - Fill CR(d_IPI) by interpolating CR vs distance (Bat-specific)
%
%   To prevent misleading extrapolation, each fill is GATED by a nearest-
%   neighbour distance criterion:
%        if min|d_query - d_support| > fill.maxGap_m  =>  set fill to NaN
%
%   i.e., we only accept a filled value if there is an actual data point
%   sufficiently close in the other dataset.
%
% -------------------------------------------------------------------------
% ANCHORING k_r FROM CR
%   k_r is estimated from near-target CR data as follows:
%     1) Take the closest anchor.closestQuantile of CR-distance points.
%     2) Let CR* be the anchor.ratePercentile percentile of CR in that subset.
%     3) Let d* be either the median or minimum distance in that subset.
%     4) Solve:
%          k_r = (c / (2*d*CR*)) - 1
%        and clamp to >= 0.
%
%   Optional bootstrap perturbs d and CR in the near-target subset to obtain
%   a 95% CI on k_r (rough sensitivity / digitisation error proxy).
%
% -------------------------------------------------------------------------
% PLOTTING CONVENTIONS
%   - Measured points: filled dots ('.')
%   - Filled (interpolated) counterpart points: open circles ('o')
%   - Fitted CR_pred(d) curve: solid line per bat
%   - Unity line in Fig 2: dashed black (d_est = d_rep)
%
% -------------------------------------------------------------------------
% USER-TUNABLE PARAMETERS (MOST IMPORTANT)
%   fill.method   : 'linear' or 'pchip' (shape-preserving cubic)
%   fill.maxGap_m : maximum allowed distance to nearest support point for fills
%   anchor.closestQuantile : fraction of closest distances used for anchoring
%   anchor.ratePercentile  : high-rate percentile used as CR*
%   doBootstrap / boot.*   : enable CI estimates for k_r (optional)
%
% -------------------------------------------------------------------------
% OUTPUTS / DIAGNOSTICS (COMMAND WINDOW)
%   For each bat:
%     - k_r estimate and CI95 (if bootstrap)
%     - median and mean ratio d_est/d_rep from measured IPI inversion
%     - number of successful fills in each direction
%
% -------------------------------------------------------------------------
% NOTES / LIMITATIONS
%   - “Fill” points are for visual UNION completeness and sanity checks; they
%     should not be interpreted as independent measurements.
%   - The nearest-neighbour gate mitigates, but does not fully eliminate,
%     interpolation bias when datasets are sparse or unevenly sampled.
%   - This script assumes monotonic distance axes are not guaranteed; data are
%     sorted and duplicate distances are removed prior to interpolation.
% -------------------------------------------------------------------------

clear; clc; close all
saveFigures = 1;

%% ---------------- USER SETTINGS ----------------
c = 343;

crFiles  = { ...
    '../data/hartley_1989/Bat1_CR_Distance.csv', ...
    '../data/hartley_1989/Bat2_CR_Distance.csv'  ...
    };

ipiFiles = { ...
    '../data/hartley_1989/Bat_1_IPI_Distance.csv', ...
    '../data/hartley_1989/Bat2_IPI_Distance.csv'   ...
    };

batNames = {'Bat 1','Bat 2'};

anchor.closestQuantile = 0.2;
anchor.ratePercentile  = 95;
anchor.useMedianD      = true;

doBootstrap = true;
boot.n = 500;
boot.sigma_d_m   = 0.02;
boot.sigma_cr_Hz = 3.0;

% Interpolation settings for "fill"
fill = struct();
fill.method = 'pchip';      % 'linear' or 'pchip'
fill.maxGap_m = 0.08;       % only fill if nearest neighbour distance <= this (avoid wild extrapolation)

% Plot ranges
dPlotMax = 5;
crPlotMax = 220;

%% ---------------- LOAD DATA ----------------
B = struct([]);
for i = 1:2
    [d_cr, cr] = localRead2col(crFiles{i});
    [d_ipi, ipi_ms] = localRead2col(ipiFiles{i});

    good1 = isfinite(d_cr) & isfinite(cr) & d_cr>0 & cr>0;
    d_cr = d_cr(good1); cr = cr(good1);

    good2 = isfinite(d_ipi) & isfinite(ipi_ms) & d_ipi>0 & ipi_ms>0;
    d_ipi = d_ipi(good2); ipi_s = ipi_ms(good2)/1000;

    B(i).name = batNames{i};
    B(i).d_cr = d_cr(:);
    B(i).cr   = cr(:);
    B(i).d_ipi = d_ipi(:);
    B(i).ipi_s = ipi_s(:);
end

colors = lines(2);

%% ---------------- BUILD INTERPOLANTS (CR(d) and IPI(d)) ----------------
for i = 1:2
    % For robustness: sort by distance and remove duplicate distances
    [dC, idx] = sort(B(i).d_cr);
    crC = B(i).cr(idx);
    [dC, ia] = unique(dC, 'stable'); crC = crC(ia);

    [dI, idx] = sort(B(i).d_ipi);
    ipiI = B(i).ipi_s(idx);
    [dI, ia] = unique(dI, 'stable'); ipiI = ipiI(ia);

    % Interpolants
    B(i).CR_of_d  = @(dq) interp1(dC, crC, dq, fill.method, NaN);
    B(i).IPI_of_d = @(dq) interp1(dI, ipiI, dq, fill.method, NaN);

    % Also store for gap checking (nearest neighbour distances)
    B(i).dC_sorted = dC;
    B(i).dI_sorted = dI;
end

%% ---------------- FIT kr FROM CR (ANCHOR) ----------------
Fits = struct([]);
for i = 1:2
    d = B(i).d_cr;
    cr = B(i).cr;

    dCut = quantile(d, anchor.closestQuantile);
    idx = d <= dCut;

    dNear = d(idx);
    crNear = cr(idx);

    crStar = prctile(crNear, anchor.ratePercentile);
    if anchor.useMedianD
        dStar = median(dNear);
    else
        dStar = min(dNear);
    end

    kr_anchor = max((c / (2*dStar*crStar)) - 1, 0);

    % Bootstrap CI (anchored estimator)
    if doBootstrap
        krBoot = nan(boot.n,1);
        for b = 1:boot.n
            dB  = dNear + boot.sigma_d_m    * randn(size(dNear));
            crB = crNear + boot.sigma_cr_Hz * randn(size(crNear));
            goodB = isfinite(dB) & isfinite(crB) & dB>0 & crB>0;
            dB = dB(goodB); crB = crB(goodB);
            if numel(dB) < 5, continue; end
            crStarB = prctile(crB, anchor.ratePercentile);
            dStarB  = anchor.useMedianD * median(dB) + (~anchor.useMedianD) * min(dB);
            krBoot(b) = max((c/(2*dStarB*crStarB))-1, 0);
        end
        ci95 = prctile(krBoot(isfinite(krBoot)), [2.5 97.5]);
    else
        ci95 = [NaN NaN];
    end

    Fits(i).name = B(i).name;
    Fits(i).kr = kr_anchor;
    Fits(i).dCut = dCut;
    Fits(i).dStar = dStar;
    Fits(i).crStar = crStar;
    Fits(i).ci95 = ci95;
end

%% ---------------- UNION-AUGMENTED ARRAYS FOR PLOTTING ----------------
U = struct([]);
for i = 1:2
    krUse = Fits(i).kr;

    % A) All CR points (measured) + filled IPI at those distances (for Fig 2)
    d_from_CR = B(i).d_cr;
    ipi_fill_at_dCR = B(i).IPI_of_d(d_from_CR);

    % gate fills by nearest-neighbour distance in the IPI domain
    nnI = localNearestDist(d_from_CR, B(i).dI_sorted);
    ipi_fill_at_dCR(nnI > fill.maxGap_m) = NaN;

    d_est_from_filledIPI = (c .* ipi_fill_at_dCR) ./ (2*(1+krUse));

    % B) All IPI points (measured) + filled CR at those distances (for Fig 1)
    d_from_IPI = B(i).d_ipi;
    cr_fill_at_dIPI = B(i).CR_of_d(d_from_IPI);

    nnC = localNearestDist(d_from_IPI, B(i).dC_sorted);
    cr_fill_at_dIPI(nnC > fill.maxGap_m) = NaN;

    % Store
    U(i).name = B(i).name;
    U(i).krUse = krUse;

    U(i).d_from_CR = d_from_CR;
    U(i).cr_meas   = B(i).cr;
    U(i).ipi_fill_at_dCR = ipi_fill_at_dCR;
    U(i).d_est_fill = d_est_from_filledIPI;

    U(i).d_from_IPI = d_from_IPI;
    U(i).ipi_meas   = B(i).ipi_s;
    U(i).cr_fill_at_dIPI = cr_fill_at_dIPI;

    % Inversion for measured IPI points (what you already do)
    U(i).d_est_meas = (c .* B(i).ipi_s) ./ (2*(1+krUse));

    % Ratios for measured IPI inversion (your key result)
    U(i).median_ratio = median(U(i).d_est_meas ./ U(i).d_from_IPI, 'omitnan');
    U(i).mean_ratio   = mean(U(i).d_est_meas ./ U(i).d_from_IPI, 'omitnan');

    % Diagnostics: how many fills succeeded?
    U(i).n_fill_IPI = sum(isfinite(ipi_fill_at_dCR));
    U(i).n_fill_CR  = sum(isfinite(cr_fill_at_dIPI));
end

%% ---------------- FIGURE 1: CR vs DISTANCE (MEASURED + FILLED) ----------------
figure('Color','w','Position',[80 80 650 320]); hold on; grid on;

for i = 1:2
    % measured CR points
    plot(U(i).cr_meas, U(i).d_from_CR, '.', 'Color', colors(i,:), 'MarkerSize', 12, ...
        'DisplayName', sprintf('%s', U(i).name));

    % filled CR values at IPI distances (open circles)
    ok = isfinite(U(i).cr_fill_at_dIPI);
    plot(U(i).cr_fill_at_dIPI(ok), U(i).d_from_IPI(ok), 'o', ...
        'Color', colors(i,:), 'MarkerSize', 4, 'LineWidth', 0.8, 'MarkerFaceColor','none', ...
        'HandleVisibility','off');
end

% overlay fitted curve
dLine = linspace(0.05, dPlotMax, 400);
for i = 1:2
    crLine = localCrPred(dLine, c, U(i).krUse);
    plot(crLine, dLine, '-', 'Color', colors(i,:), 'LineWidth', 2.0, ...
        'DisplayName', sprintf('%s fit: $k_r$=%.2f', U(i).name, U(i).krUse));
end

xlabel('Call rate (Hz)','Interpreter','latex');
ylabel('Reported target distance (m)','Interpreter','latex');
title('Hartley et al. (1989) reanalysis: $k_r$ anchored from $C_r$', 'Interpreter','latex');
set(gca,'YLim',[0 dPlotMax], 'XLim',[0 crPlotMax], 'TickLabelInterpreter','latex', 'FontSize',11);
legend('Location','northeast','Interpreter','latex');
formatLatex(gca)
if saveFigures
    saveFigure(gcf, '../fig/', 'hartely_cr_dist_responsivity_UNION')
end

%% ---------------- FIGURE 2: IPI-INVERTED DISTANCE vs REPORTED (MEASURED + FILLED) ----------------
figure('Color','w','Position',[120 120 650 320]); hold on; grid on;

for i = 1:2
    % measured IPI inversion points
    plot(U(i).d_from_IPI, U(i).d_est_meas, '.', 'Color', colors(i,:), 'MarkerSize', 12, ...
        'DisplayName', sprintf('%s', U(i).name));

    % filled IPI-at-CR-distances inversion points (open circles)
    ok = isfinite(U(i).d_est_fill);
    plot(U(i).d_from_CR(ok), U(i).d_est_fill(ok), 'o', ...
        'Color', colors(i,:), 'MarkerSize', 4, 'LineWidth', 0.8, 'MarkerFaceColor','none', ...
        'HandleVisibility','off');
end

lims = [0 dPlotMax];
plot(lims, lims, 'k--', 'LineWidth', 1.5, 'DisplayName','unity $d_{est}=d_{rep}$');

xlabel('Reported target distance (m)','Interpreter','latex');
ylabel('Cue distance implied by IPI (m)','Interpreter','latex');
title('Inverting IPI with $C_r$-anchored $k_r$', 'Interpreter','latex');
set(gca,'XLim',lims,'YLim',lims,'TickLabelInterpreter','latex','FontSize',11);
legend('Location','northwest','Interpreter','latex');
formatLatex(gca)
if saveFigures
    saveFigure(gcf, '../fig/', 'hartely_IPI_dist_responsivity_UNION')
end

%% ---------------- PRINT SUMMARY ----------------
fprintf('\n=== Hartley Fig9 (UNION / fill-by-distance): c=%g m/s | fill.maxGap=%.3f m ===\n', ...
    c, fill.maxGap_m);

for i = 1:2
    fprintf('%s:\n', U(i).name);
    fprintf('  kr_anchor = %.3f | CI95=[%.3f, %.3f]\n', Fits(i).kr, Fits(i).ci95(1), Fits(i).ci95(2));
    fprintf('  IPI inversion (measured IPI points): median(d_est/d_rep)=%.3f | mean=%.3f\n', ...
        U(i).median_ratio, U(i).mean_ratio);
    fprintf('  fills: IPI@CR-dist succeeded for %d/%d CR points\n', U(i).n_fill_IPI, numel(U(i).d_from_CR));
    fprintf('        CR@IPI-dist succeeded for %d/%d IPI points\n', U(i).n_fill_CR,  numel(U(i).d_from_IPI));
end

%% ---------------- LOCAL HELPERS ----------------
function crPred = localCrPred(d, c, kr)
    crPred = c ./ (2*(1+kr).*d);
end

function [x, y] = localRead2col(fname)
    if ~isfile(fname)
        error('File not found: %s', fname);
    end
    M = readmatrix(fname);
    if size(M,2) < 2
        T = readtable(fname);
        M = T{:,:};
    end
    x = M(:,1);
    y = M(:,2);
end

function nn = localNearestDist(query_d, support_d_sorted)
    % nearest-neighbour distance from each query point to a sorted support array
    query_d = query_d(:);
    support_d_sorted = support_d_sorted(:);
    nn = nan(size(query_d));

    for k = 1:numel(query_d)
        [m, ~] = min(abs(support_d_sorted - query_d(k)));
        nn(k) = m;
    end
end