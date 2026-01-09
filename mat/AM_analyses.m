%% Batch AM & Ripple analysis for all WAV files (per channel)
% - Scans rootDir recursively for CALL_*.wav
% - Per file/channel:
%     • Envelope via Hilbert + low-pass
%     • Duration via percentdur (default 10–90% energy)
%     • fAM from envelope (Welch; autocorr fallback)
%     • Echo-delay metrics (xcorr/cepstrum/comb)  →  fRipple = 1000/τ_cons(ms)
% - Writes CSV and makes summary plots (incl. Duration vs fAM and Duration vs fRipple)

clear; clc;

%% ===== USER SETTINGS =====
rootDir     = '../data/'; % large dataset is not provided due to upload limits.
outCSV      = fullfile('results', 'am_analysis_results.csv');
makeSummary = true;                      % quick plots at the end
recurse     = true;                      % scan subfolders
fileGlob    = 'CALL_*.wav';              % only files starting with "CALL_"

% Envelope & analysis params
envLP_Hz    = 5e3;                       % low-pass for envelope (Hz)
durPct      = [0.10 0.90];               % duration window from cum. energy
fAM_band    = [300 20000];               % Hz, search band for AM peak
minCallMs   = 0.1;                       % ignore < this after windowing (ms)
useEnvForDur= false;                     % derive duration window from envelope energy

% Welch spectrum params (envelope)
welchNFFT   = 4096;
welchWlen   = 512;
welchOv     = 0.5;

% Quality filter
peakThreshold = 0.0035;                  % reject weak calls

%% ===== FIND FILES =====
if recurse
    files = dir(fullfile(rootDir, '**', fileGlob));   % R2016b+
else
    files = dir(fullfile(rootDir, fileGlob));
end
files = files(~[files.isdir]);
if isempty(files)
    warning('No WAV files found under: %s', rootDir);
end

%% ===== RESULT TABLE INIT (now includes fRipple) =====
Results = table('Size',[0 25], ...
    'VariableTypes', {'string','string','double','double','double', ...
                      'double','double','double','double','double', ...
                      'double','double','double','double','string', ...
                      'double','double','double','double','double','double','double', ...
                      'double','double','double'}, ...
    'VariableNames', {'File','Folder','Fs','NumChannels','Channel', ...
                      'Dur_ms','tL_ms','tR_ms','fAM_Hz','fAM_kHz', ...
                      'PeakProminence','EnvLP_Hz','Presence_RMS','Fmax_kHz','AM_Method', ...
                      'Tau_xcorr_ms','Tau_cep_ms','Tau_ripple_ms','Tau_cons_ms', ...
                      'OverlapFlag','nCorrPeaks','DurInfl_ms', ...
                      'fRipple_Hz','fRipple_kHz','DurCorr_ms'});

%% ===== PROCESS LOOP =====
nFiles = numel(files);
tStart = tic;

for i = 1:nFiles
    fpath = fullfile(files(i).folder, files(i).name);

    % progress tracker
    if mod(i,10)==1 || i==nFiles
        elapsed = toc(tStart);
        avgTime = elapsed / max(i,1);
        eta = avgTime * (nFiles - i);
        fprintf('Processing %d/%d (%.1f%%) | Elapsed: %.1fs | ETA: %.1fs\n', ...
            i, nFiles, 100*i/nFiles, elapsed, max(eta,0));
    end

    % read
    try
        [x, fs] = audioread(fpath);
    catch ME
        warning('Failed to read %s: %s', fpath, ME.message);
        continue;
    end
    if isempty(x)
        warning('Empty audio in %s (skipping).', fpath);
        continue;
    end

    % quality filter
    if max(abs(x),[],'all') < peakThreshold
        fprintf('Skipping %s (peak below threshold)\n', files(i).name);
        continue;
    end

    % basics
    nCh  = size(x,2);
    N    = size(x,1);
    t_ms = (0:N-1)/fs*1000;

    % envelope LP
    Wc = min(envLP_Hz, 0.99*(fs/2)) / (fs/2);
    if Wc <= 0
        warning('envLP_Hz too low relative to fs in %s; skipping.', fpath);
        continue;
    end
    [bLP, aLP] = butter(4, Wc, 'low');

    % per-channel analysis
    for ch = 1:nCh
        sig = x(:,ch);

        % Envelope
        sigEnv = abs(hilbert(sig));
        sigEnv = filtfilt(bLP, aLP, sigEnv);

        % Duration from your percentdur (on waveform)
        [durIdx, pres, fmax_kHz] = percentdur(sig, durPct);
        iL = max(1, min(durIdx(1), N));
        iR = max(1, min(durIdx(2), N));
        leftTime  = t_ms(iL);
        rightTime = t_ms(iR);
        dur_ms    = rightTime - leftTime;

        % AM freq from envelope within window
        segEnv = sigEnv(iL:iR);
        if ~any(segEnv)
            fAM = NaN; prom = NaN; method = "none";
        else
            segN = segEnv ./ max(max(segEnv), eps);
            [fAM, prom, method] = estimate_fAM_from_envelope( ...
                segN, fs, fAM_band, welchNFFT, welchWlen, welchOv);
        end

        % Echo metrics → ripple (Hz)
        rawSeg = sig(iL:iR);
        echoOpts = struct( ...
            'band_Hz',    [20e3 90e3], ...
            'tau_ms_rng', [0.2 2.0], ...
            'wlen_ms',    8, ...
            'xc_minprom', 0.10, ...
            'comb_minR',  0.30, ...
            'agree_ms',   0.15);
        M = detect_overlap_metrics(rawSeg, fs, echoOpts);

        % Duration inflation heuristic
        durInfl_ms   = M.is_overlapped * M.consensus_ms;

        % fRipple from consensus tau (guard NaN/Inf)
        if ~isnan(M.consensus_ms) && M.consensus_ms > 0
            fRipple_Hz  = 1000 ./ M.consensus_ms;
        else
            fRipple_Hz  = NaN;
        end

                % Duration from your percentdur (on waveform)
        [durIdx, pres, fmax_kHz] = percentdur(sig, durPct);
        iL = max(1, min(durIdx(1), N));
        iR = max(1, min(durIdx(2), N));
        leftTime  = t_ms(iL);
        rightTime = t_ms(iR);
        dur_ms    = rightTime - leftTime;

        % AM freq from envelope within window
        segEnv = sigEnv(iL:iR);
        if ~any(segEnv)
            fAM = NaN; prom = NaN; method = "none";
        else
            segN = sigEnv ./ max(max(sigEnv), eps);
            [fAM, prom, method] = estimate_fAM_from_envelope( ...
                segN, fs, fAM_band, welchNFFT, welchWlen, welchOv);
        end

        % Echo metrics → ripple (Hz)
        rawSeg = sig(iL:iR);
        M = detect_overlap_metrics(rawSeg, fs, echoOpts);

        % Duration inflation heuristic
        durInfl_ms   = M.is_overlapped * M.consensus_ms;

        % Corrected duration (clip at >=0)
        durCorr_ms = dur_ms;
        if ~isnan(durInfl_ms)
            durCorr_ms = max(0, dur_ms - durInfl_ms);
        end

        % fRipple from consensus tau
        if ~isnan(M.consensus_ms) && M.consensus_ms > 0
            fRipple_Hz  = 1000 ./ M.consensus_ms;
        else
            fRipple_Hz  = NaN;
        end

        % Append row
        Results = [Results; {
            string(files(i).name), string(files(i).folder), fs, nCh, ch, ...
            dur_ms, leftTime, rightTime, fAM, fAM/1000, prom, envLP_Hz, ...
            pres, fmax_kHz, string(method), ...
            M.tau_xcorr_ms, M.tau_cep_ms, M.tau_ripple_ms, M.consensus_ms, ...
            M.is_overlapped, double(M.n_corr_peaks), durInfl_ms, ...
            fRipple_Hz, fRipple_Hz/1000, durCorr_ms}]; %#ok<AGROW>
    end
end

%% ===== SAVE RESULTS =====
if ~isempty(Results)
    writetable(Results, outCSV);
    fprintf('Saved results to: %s\n', outCSV);
else
    warning('No results to save.');
end

%% ===== SUMMARY PLOTS (per-channel histograms + correlations) =====
if makeSummary && ~isempty(Results)
    figure('Color','w','Name','AM & Ripple Summary', 'position', [100 100 1400 600]);
    tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

    % === Duration histogram (per channel) ===
    nexttile; hold on;
    channels = unique(Results.Channel);
    cmap = lines(numel(channels));
    for ci = 1:numel(channels)
        v = Results.Dur_ms(Results.Channel==channels(ci));
        h = histogram(v,40,'Normalization','count','DisplayStyle','bar');
        h.FaceColor = cmap(ci,:);
        h.FaceAlpha = 0.35;
        h.EdgeAlpha = 0.5;
    end
    grid on; xlabel('Duration (ms)'); ylabel('Count');
    title('Call duration');
    % legend("Ch "+string(channels),'Location','best'); legend boxoff;
    formatLatex(gca);

    % === fAM histogram (per channel) ===
    nexttile; hold on;
    for ci = 1:numel(channels)
        v = Results.fAM_Hz(Results.Channel==channels(ci));
        h = histogram(v,40,'Normalization','count','DisplayStyle','bar');
        h.FaceColor = cmap(ci,:);
        h.FaceAlpha = 0.35;
        h.EdgeAlpha = 0.5;
    end
    grid on; xlabel('$f_{AM}$ (Hz)'); ylabel('Count');
    title('AM frequency');
    % legend("Ch "+string(channels),'Location','best'); legend boxoff;
    formatLatex(gca);

    % === fRipple histogram (per channel) ===
    nexttile; hold on;
    for ci = 1:numel(channels)
        v = Results.fRipple_Hz(Results.Channel==channels(ci));
        h = histogram(v,40,'Normalization','count','DisplayStyle','bar');
        h.FaceColor = cmap(ci,:);
        h.FaceAlpha = 0.35;
        h.EdgeAlpha = 0.5;
    end
    grid on; xlabel('$f_{ripple}$ (Hz)'); ylabel('Count');
    title('Ripple frequency');
    legend("Ch "+string(channels),'Location','best', 'Interpreter','latex'); legend boxoff;
    formatLatex(gca);

    % === Duration vs fAM (+ correlation text) ===
    nexttile;
    v = ~isnan(Results.Dur_ms) & ~isnan(Results.fAM_Hz);
    scatter(Results.Dur_ms(v), Results.fAM_Hz(v), 12, 'filled', 'MarkerFaceAlpha',0.35);
    box on; xlabel('Duration (ms)'); ylabel('$f_{AM}$ (Hz)'); title('$f_{AM}$ vs Duration');
    if any(v)
        [rhoA,pA] = corr(Results.Dur_ms(v), Results.fAM_Hz(v), 'Type','Spearman');
        [rA,pA2]  = corr(Results.Dur_ms(v), Results.fAM_Hz(v), 'Type','Pearson');
        xl = xlim; yl = ylim;
        text(xl(1)+0.01*range(xl), yl(2)-0.01*range(yl), ...
            sprintf('$\\rho=%.2f (p=%.1e)$ \\ $r=%.2f (p=%.1e)$', rhoA,pA,rA,pA2), ...
            'BackgroundColor','w','EdgeColor','k','FontSize',11,'interpreter','latex');
    end
    formatLatex(gca);

    % === Duration vs fRipple (+ correlation text) ===
    nexttile;
    v = ~isnan(Results.Dur_ms) & ~isnan(Results.fRipple_Hz);
    scatter(Results.Dur_ms(v), Results.fRipple_Hz(v), 12, 'filled', 'MarkerFaceAlpha',0.35);
    grid on; box on; xlabel('Duration (ms)'); ylabel('$f_{ripple}$ (Hz)'); title('$f_{ripple}$ vs Duration');
    if any(v)
        [rhoR,pR] = corr(Results.Dur_ms(v), Results.fRipple_Hz(v), 'Type','Spearman');
        [rR,pR2]  = corr(Results.Dur_ms(v), Results.fRipple_Hz(v), 'Type','Pearson');
        xl = xlim; yl = ylim;
        text(xl(1)+0.1*range(xl), yl(2)-0.95*range(yl), ...
            sprintf('$\\rho=%.2f (p=%.1e) \\ r=%.2f (p=%.1e)$', rhoR,pR,rR,pR2), ...
            'BackgroundColor','w','EdgeColor','k','FontSize',11,'interpreter','latex');
    end
    formatLatex(gca);

    % === Per-channel fAM boxplot (as-is) ===
    nexttile;
    boxchart(categorical(Results.Channel), Results.fAM_Hz);
    xlabel('Channel'); ylabel('$f_{AM}$ (Hz)'); title('Per-channel $f_{AM}$');
    formatLatex(gca);
end

%%
meanDur = mean(Results.Dur_ms,'omitnan');
stdDur  = std(Results.Dur_ms,'omitnan');

explained_var = eta2 * stdDur^2;   % variance explained by channel
explained_sd  = sqrt(explained_var);

fprintf('Mean duration = %.3f ms\n', meanDur);
fprintf('Total SD = %.3f ms\n', stdDur);
fprintf('Explained SD (channel effect) ≈ %.3f ms\n', explained_sd);
fprintf('Explained SD as %% of mean ≈ %.1f%%\n', 100*explained_sd/meanDur);


%% Save Figure
fig_path = '../manuscript/fig';
base_name = 'modulation_analyses';
saveFigure(gcf, fig_path, base_name);
%% ===== Channelwise statistical analysis =====
channels = unique(Results.Channel);

% --- Parameters to test ---
params = {
    'Dur_ms',      'Duration (ms)';
    'DurCorr_ms',  'Corrected Duration (ms)';
    'fAM_Hz',      'AM frequency (Hz)';
    'fRipple_Hz',  'Ripple frequency (Hz)';
};

StatsTable = table('Size',[0 6], ...
    'VariableTypes',{'string','double','double','double','double','double'}, ...
    'VariableNames',{'Parameter','pKW','H','N','eta2','eps2'});

for pi = 1:size(params,1)
    paramName  = params{pi,1};
    paramLabel = params{pi,2};

    % Extract valid
    vals = Results.(paramName);
    idx  = ~isnan(vals);
    vals = vals(idx);
    ch_labels = Results.Channel(idx);

    % --- Kruskal–Wallis test ---
    [pKW, tbl, stats] = kruskalwallis(vals, ch_labels, 'off');

    % Extract effect size
    H = tbl{2,5};    
    N = sum(cell2mat(tbl(2:end-1,3))); 
    k = numel(unique(ch_labels));

    eta2 = (H - (k-1)) / (N - 1);
    eps2 = (H - (k-1)) / (N + k - 1);

    % Store
    StatsTable = [StatsTable; {paramLabel, pKW, H, N, eta2, eps2}];

    % Print to console
    fprintf('\n--- %s ---\n', paramLabel);
    fprintf('Kruskal–Wallis p = %.4g, H = %.3f, N = %d\n', pKW, H, N);
    fprintf('Effect size: η² = %.4f, ε² = %.4f\n', eta2, eps2);

    % --- Post-hoc if significant ---
    if pKW < 0.05
        figure('Color','w');
        multcompare(stats, 'CType','dunn-sidak');
        title(sprintf('Post-hoc channel comparisons (%s)', paramLabel));
    end
end

% ===== Correlations =====
% Duration vs fAM
valid = ~isnan(Results.Dur_ms) & ~isnan(Results.fAM_Hz);
[rho1, p1] = corr(Results.Dur_ms(valid), Results.fAM_Hz(valid), 'Type','Spearman');
[r1, p1_lin] = corr(Results.Dur_ms(valid), Results.fAM_Hz(valid), 'Type','Pearson');

% Duration vs fRipple
valid = ~isnan(Results.Dur_ms) & ~isnan(Results.fRipple_Hz);
[rho2, p2] = corr(Results.Dur_ms(valid), Results.fRipple_Hz(valid), 'Type','Spearman');
[r2, p2_lin] = corr(Results.Dur_ms(valid), Results.fRipple_Hz(valid), 'Type','Pearson');

% Append correlation results to table
CorrTable = table(...
    ["Duration vs fAM"; "Duration vs fRipple"], ...
    [rho1; rho2], [p1; p2], [r1; r2], [p1_lin; p2_lin], ...
    'VariableNames',{'Comparison','Spearman_rho','Spearman_p','Pearson_r','Pearson_p'});

% ===== Final Output =====
disp('===== Channelwise Kruskal–Wallis Results =====');
disp(StatsTable);

disp('===== Correlation Results =====');
disp(CorrTable);
%% === Turn off caffeinate at the end ===
% if ~isempty(caffeinatePID) && caffeinatePID > 0
%     system(sprintf('kill %d', caffeinatePID));
%     fprintf('caffeinate process %d killed.\n', caffeinatePID);
% end
%% =================== HELPERS ===================

function [iL, iR] = duration_idx_energy(sig, pct)
% Indices enclosing [pct(1) .. pct(2)] of cumulative energy
    s = sig(:);
    E = cumsum(s.^2);
    Et = E(end);
    if Et <= eps
        iL = 1; iR = min(2, numel(s)); return;
    end
    E = E / Et;

    % gentle smoothing to avoid noisy bounds
    sm = max(5, 2*floor((numel(s)/10)/2)+1); % odd, ~10% len
    E = conv(E, ones(sm,1)/sm, 'same');

    p1 = max(0, min(pct(1), pct(2)));
    p2 = min(1, max(pct(1), pct(2)));

    iL = find(E >= p1, 1, 'first');
    iR = find(E >= p2, 1, 'first');
    if isempty(iL), iL = 1; end
    if isempty(iR), iR = numel(s); end
end

function [fAM, peakProm, method] = estimate_fAM_from_envelope(envNorm, fs, bandHz, nfft, wlen, ovFrac)
% estimate_fAM_from_envelope:
%   Primary: Welch spectrum peak in band
%   Fallback: autocorr peak (positive lag)
    method   = 'welch';
    peakProm = NaN;
    fAM      = NaN;

    % Detrend & zero-mean
    seg = detrend(envNorm(:), 'linear');
    seg = seg - mean(seg);

    % Welch spectrum
    if nargin < 4 || isempty(nfft), nfft = 4096; end
    if nargin < 5 || isempty(wlen), wlen = 512;  end
    if nargin < 6 || isempty(ovFrac), ovFrac = 0.5; end
    nover = round(ovFrac*wlen);

    [P, f] = pwelch(seg, hann(wlen), nover, nfft, fs);
    mask   = (f >= bandHz(1)) & (f <= bandHz(2));

    if any(mask)
        [pks, locs] = findpeaks(P(mask), f(mask), 'SortStr','descend');
        if ~isempty(pks)
            fAM      = locs(1);
            peakProm = pks(1) / max(P(mask)+eps);
        end
    end

    % Fallback: autocorrelation if needed
    if isnan(fAM)
        method = 'autocorr';
        [r,lags] = xcorr(seg, 'coeff');
        mid = find(lags==0,1);
        lagsP = lags(mid+1:end);
        rP    = r(mid+1:end);

        % guard tiny lags (avoid spuriously high f)
        guard = max(1, round(0.2e-3*fs)); % 0.2 ms
        if numel(rP) > guard
            [~, idx] = max(rP(guard:end));
            lagSamp = lagsP(guard-1+idx);
            if lagSamp > 0
                fAM = fs / lagSamp;
                peakProm = rP(guard-1+idx);
            end
        end
    end
end

% function [duration,pres,fmax] = percentdur(seg,frac)
% % duration: [iL iR] sample indices covering frac energy (e.g., [0.10 0.90])
% % pres    : RMS within [iL iR]
% % fmax    : dominant bin (kHz) from an NFFT=192 FFT → 1 kHz/bin at fs=192 kHz
% 
%     sm = 2*floor((length(seg)/10)/2)+1;     % odd length ~10% of signal
%     csig = cumsum(seg.^2);
%     csig(1) = 0;
%     csig = csig/max(csig);
%     csig = smooth(smooth(csig,sm),sm);
% 
%     int = find(csig>frac(1) & csig<frac(2));
%     if isempty(int)
%         duration = [1 2];
%         pres = 0;
%         fmax = NaN;
%     else
%         duration = [int(1) int(end)];
%         pres = rms(seg(int));
%         mag = 20*log10(abs(fft(seg(int(1):int(end)), 192))); % 1 kHz/bin @ fs=192 kHz
%         [~, fmax] = max(mag(1:length(mag)/2));
%         fmax = fmax - 1;   % convert bin index → kHz (bin 1 = 0 kHz)
%     end
% end

function M = detect_overlap_metrics(sig, fs, opts)
% DETECT_OVERLAP_METRICS  Multipath/overlap detector for bat calls
% Returns echo-delay estimates from:
%   - time-domain autocorr/xcorr (tau_xcorr)
%   - real-cepstrum (tau_cep)
%   - spectral-comb spacing (tau_ripple)
% and a consensus flag.
%
% Inputs
%   sig  : column vector (double)
%   fs   : sample rate (Hz)
%   opts : struct with fields (all optional):
%          .band_Hz     = [20e3 90e3]    % analysis band for spectrum tests
%          .tau_ms_rng  = [0.2 2.0]      % plausible echo lags (ms)
%          .wlen_ms     = 8              % analysis window length (ms) if trimming around peak
%          .xc_minprom  = 0.1            % min normalized peak prominence (0-1) in xcorr
%          .comb_minR   = 0.3            % min autocorr of spectrum to accept a spacing
%          .agree_ms    = 0.15           % max pairwise disagreement (ms) to call consensus
%
% Output struct M:
%   .tau_xcorr_ms, .tau_cep_ms, .tau_ripple_ms
%   .n_corr_peaks    (# strong xcorr peaks in positive lags within tau range)
%   .consensus_ms    (median of available taus)
%   .is_overlapped   (logical)
%   .notes           (string)

% --- defaults
if nargin < 3, opts = struct; end
fld = @(f,d) (isfield(opts,f) && ~isempty(opts.(f))) * opts.(f) + (~(isfield(opts,f)&&~isempty(opts.(f))))*d;
band_Hz    = fld('band_Hz',    [20e3 90e3]);
tau_ms_rng = fld('tau_ms_rng', [0.2 2.0]);
wlen_ms    = fld('wlen_ms',    8);
xc_minprom = fld('xc_minprom', 0.10);
comb_minR  = fld('comb_minR',  0.30);
agree_ms   = fld('agree_ms',   0.15);

sig = sig(:);
sig(~isfinite(sig)) = 0;
if ~any(sig), M = emptyOut(); M.notes="empty sig"; return; end

% ---- optional trim around global peak to stabilize analysis
N  = numel(sig);
pk = max(abs(sig));
if pk > 0
    [~, imax] = max(abs(sig));
    half = round((wlen_ms/1000)*fs/2);
    i1 = max(1, imax-half); i2 = min(N, imax+half);
    s = sig(i1:i2);
else
    s = sig;
end

% Detrend
s = detrend(s, 'linear'); s = s - mean(s);

% ========== 1) xcorr in time ==========
[r,lags] = xcorr(s, 'coeff');
mid   = find(lags==0,1);
lagsP = lags(mid+1:end) / fs;               % seconds, positive lags
rP    = r(mid+1:end);

% search only in plausible tau window
tau_lo = tau_ms_rng(1)/1000; tau_hi = tau_ms_rng(2)/1000;
mask   = (lagsP >= tau_lo) & (lagsP <= tau_hi);

tau_xcorr_ms = NaN; n_corr_peaks = 0;
if any(mask)
    rp  = rP(mask);
    lg  = lagsP(mask);
    [pks,locs] = findpeaks(rp, 'MinPeakProminence', xc_minprom);
    n_corr_peaks = numel(pks);
    if ~isempty(pks)
        [~,ix] = max(pks);
        tau_xcorr_ms = lg(locs(ix))*1000;   % ms
    end
end

% ========== 2) real cepstrum ==========
% real cepstrum of magnitude spectrum; peak at quefrency ~ tau
w  = hann(numel(s));
S  = abs(fft(s.*w));
S  = S + eps;
C  = real(ifft(log(S)));
q  = (0:numel(C)-1).'/fs;                  % seconds (quefrency)
% restrict to plausible region
cmask = (q >= tau_lo) & (q <= tau_hi);
tau_cep_ms = NaN;
if any(cmask)
    [cpk,cloc] = max(C(cmask));
    if ~isempty(cpk)
        qv = q(cmask); tau_cep_ms = qv(cloc)*1000;
    end
end

% ========== 3) spectral comb spacing ==========
% band-limit, then estimate ripple spacing via spectrum autocorr
[bBP,aBP] = butter(4, band_Hz/(fs/2), 'bandpass');
sb = filtfilt(bBP,aBP, s);
NFFT = 8192;
Sf   = abs(fft(sb.*hann(numel(sb)), NFFT));
f    = (0:NFFT-1)'/NFFT*fs;
% keep analysis band
bmask = (f>=band_Hz(1)) & (f<=band_Hz(2));
SfB   = Sf(bmask);
SfB   = SfB / max(SfB + eps);
% autocorr of magnitude spectrum to find periodic spacing
[Rf, lagsF] = xcorr(SfB - mean(SfB), 'coeff');
midF = find(lagsF==0,1);
RfP  = Rf(midF+1:end);
% convert spectral-lag to frequency spacing (bins → Hz)
df   = fs/NFFT;
freqSpacings_Hz = (1:numel(RfP)).' * df;
% plausible spacing from tau window
dflow = 1/(tau_hi); dfhi = 1/(tau_lo);
maskS = (freqSpacings_Hz >= dflow) & (freqSpacings_Hz <= dfhi);
tau_ripple_ms = NaN;
if any(maskS)
    [pkS,ixS] = max(RfP(maskS));
    if ~isempty(pkS) && pkS >= comb_minR
        df_star = freqSpacings_Hz(maskS); df_star = df_star(ixS);
        tau_ripple_ms = 1000/df_star;
    end
end

% ========== consensus & flags ==========
taus = [tau_xcorr_ms, tau_cep_ms, tau_ripple_ms];
taus = taus(~isnan(taus));
consensus_ms = NaN; is_overlap = false; note = "";

if numel(taus) >= 2
    med = median(taus);
    if all(abs(taus - med) <= agree_ms)
        consensus_ms = med; is_overlap = true;
        note = "agree(xcorr,cep/comb)";
    else
        % weaker: accept best pair if any two within agree_ms
        D = abs(taus(:)-taus(:).'); D(eye(size(D))==1) = Inf;
        if any(D(:) <= agree_ms)
            consensus_ms = med; is_overlap = true;
            note = "pairwise agreement";
        else
            note = "no agreement";
        end
    end
elseif numel(taus)==1
    consensus_ms = taus(1); is_overlap = true;
    note = "single method";
else
    note = "no estimate";
end

% package
M.tau_xcorr_ms  = tau_xcorr_ms;
M.tau_cep_ms    = tau_cep_ms;
M.tau_ripple_ms = tau_ripple_ms;
M.n_corr_peaks  = n_corr_peaks;
M.consensus_ms  = consensus_ms;
M.is_overlapped = is_overlap;
M.notes         = note;
M.opts_used     = struct('band_Hz',band_Hz,'tau_ms_rng',tau_ms_rng,...
                         'wlen_ms',wlen_ms,'xc_minprom',xc_minprom,...
                         'comb_minR',comb_minR,'agree_ms',agree_ms);
end

function M = emptyOut()
M = struct('tau_xcorr_ms',NaN,'tau_cep_ms',NaN,'tau_ripple_ms',NaN,...
           'n_corr_peaks',0,'consensus_ms',NaN,'is_overlapped',false,...
           'notes',"");
end