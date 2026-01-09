%% spatial_invariance_multichannel_call_characterisation.m
% -------------------------------------------------------------------------
% PURPOSE
%   Screen a large multichannel (4-ch) bat-call dataset for high-amplitude
%   segments, align channels in time, and generate a per-file diagnostic
%   “row” of plots that characterise cross-channel similarity:
%
%     (1) Normalised power spectra (pwelch) for all 4 channels
%     (2) Spectrogram (img_spec) of aligned channel 1
%     (3) Time-domain waveforms of all aligned channels + duration markers
%     (4) Normalised AM envelopes of all aligned channels + similarity metrics
%
%   The intent is to test spatial invariance of spectral ripples / envelope
%   structure across microphones (i.e., receiver-independent patterns).
%
% -------------------------------------------------------------------------
% INPUT DATA (NOT IN REPO)
%   rootDir points to a large external dataset omitted from version control:
%       rootDir = '/data/'
%
%   Files are discovered with:
%       BAT*_Analysis/segment_*/Call_*.wav
%
%   Expected audio format:
%       - .wav with 4 channels (Nx4). Files that are not 4-ch are ignored.
%
% -------------------------------------------------------------------------
% HIGH-AMPLITUDE FILE SELECTION (FIRST PASS)
%   The script first scans all candidate files and computes a global peak:
%       peakVal = max(abs(data), [], 'all')
%
%   Files are retained if:
%       peakVal > peakThreshold
%
%   Retained files are sorted by peak descending and the top N are processed:
%       maxFiles = 50
%
%   This provides a consistent “worst-case / strongest-signal” subset for
%   visual inspection and metric computation.
%
% -------------------------------------------------------------------------
% PREPROCESSING PER FILE (SECOND PASS)
%   Given a selected 4-ch file:
%
%   1) ALIGN channels 2–4 to channel 1 (reference) by maximising the
%      normalised cross-correlation:
%          [xc,lags] = xcorr(ch, ref, 'coeff');
%          lag = lags(argmax(xc));
%          ch_aligned = circshift(ch, -lag);
%
%      Note: circshift wraps around; if wrap artefacts matter, replace with
%      a zero-padded shift.
%
%   2) CENTRE each channel’s maximum absolute peak at the waveform midpoint:
%          midpoint = round(N/2);
%          shift = midpoint - peakIdx;
%          ch_centered = circshift(ch, shift);
%
%      This makes waveform and envelope comparisons visually consistent.
%
% -------------------------------------------------------------------------
% PLOT LAYOUT
%   A single figure is created with maxFiles rows and 4 columns:
%
%     Column 1: Normalised spectrum (pwelch) for all 4 channels
%     Column 2: Spectrogram (img_spec) of aligned channel 1
%     Column 3: Aligned waveforms for all 4 channels + duration markers
%     Column 4: Normalised AM envelopes for all 4 channels + summary metrics
%
%   Each row corresponds to one file from the “top-by-peak” list.
%
% -------------------------------------------------------------------------
% COLUMN 1: NORMALISED SPECTRA
%   - Uses pwelch with:
%       nfft = 1024
%       win = hann(nfft)
%       noverlap = 0.5*nfft
%
%   - Spectra are normalised per-channel to unit peak:
%       Pxx = Pxx ./ max(Pxx)
%
%   - Plotted as 10*log10(Pxx) in dB with x-axis in kHz.
%
%   Plot limits:
%       xlim([20 80]) kHz
%       ylim([-40 10]) dB (normalised)
%
% -------------------------------------------------------------------------
% COLUMN 2: SPECTROGRAM (CHANNEL 1)
%   - Uses img_spec(alignedData(:,1), 128, 127, 4096, fs, 100)
%   - Frequency tick labels converted to kHz for display.
%
%   A try/catch guards against failures in img_spec for any file.
%
% -------------------------------------------------------------------------
% COLUMN 3: ALIGNED WAVEFORMS + DURATION ESTIMATES
%   - Plots time-domain waveforms for all channels over a restricted window:
%       t in ms; xlim([5 13])
%
%   - Computes energy-based duration bounds per channel using percentdur:
%       [durationIdx,~,~] = percentdur(sig, [0.1 0.90])
%
%     This returns indices bracketing the central 10–90% cumulative energy.
%     Vertical dashed lines mark left/right times for each channel.
%
%   - Stores per-channel duration (ms) in a durationTable:
%       Dur1..Dur4 = rightTime - leftTime
%
% -------------------------------------------------------------------------
% COLUMN 4: AM ENVELOPES + CROSS-CHANNEL COMPARISONS
%   Envelope extraction:
%     env = abs(hilbert(sig))
%     env is low-pass filtered to 5 kHz:
%         butter(4, 5000/(fs/2), 'low') + filtfilt
%
%   Normalisation:
%     envNorm = env / max(env)
%
%   Pairwise comparisons for all channel pairs (6 pairs total):
%     - Pearson correlation (r)
%     - RMSE
%     - Max normalised cross-correlation and its lag (ms)
%     - DTW distance (dtw), if available (else NaN)
%
%   Results are stored in amSummaryTable with one row per file×pair.
%
% -------------------------------------------------------------------------
% CONFIGURATION PARAMETERS (TOP OF SCRIPT)
%   rootDir        : dataset root (external)
%   maxFiles       : number of top-peak files to process (default 50)
%   peakThreshold  : screening threshold on max(abs(wav)) (default 0.004)
%   chColors       : channel colour codes for plotting
%   fig_path       : output folder for saving (saveFigure is commented)
%
% -------------------------------------------------------------------------
% OUTPUTS
%   - Figure: “Aligned Spectra, Spectrograms & Waveforms”
%     Overall title:
%       'Cross-Channel Signal Characterisation - Spatial Invariance'
%
%   - durationTable (in workspace):
%       per-file duration estimates for channels 1–4 (ms)
%
%   - amSummaryTable (in workspace; optional CSV export commented):
%       per-file, per-channel-pair envelope similarity metrics
%
% -------------------------------------------------------------------------
% DEPENDENCIES
%   - Signal Processing Toolbox (pwelch, butter, filtfilt, hilbert, xcorr)
%   - dtw (Signal Processing / DSP depending on MATLAB version)
%   - percentdur 
%   - img_spec (custom helper)
%   - formatLatex (custom helper)
%   - saveFigure (custom helper; currently commented)
%
% -------------------------------------------------------------------------
rootDir = '/data/'; % large dataset, not provided in the repo due to size limits
files = dir(fullfile(rootDir, 'BAT*_Analysis', 'segment_*', 'Call_*.wav'));
fig_path = '../fig';
chColors = {'r', 'g', 'b', 'k'};
maxFiles = 50;
peakThreshold = 0.004;

% ===== First pass: collect peak values =====
filePeaks = []; % struct array: file index + peak

for i = 1:length(files)
    filePath = fullfile(files(i).folder, files(i).name);
    [data, ~] = audioread(filePath);
    if size(data,2) == 4
        peakVal = max(abs(data), [], 'all');
        if peakVal > peakThreshold
            filePeaks(end+1).index = i; %#ok<*SAGROW>
            filePeaks(end).peak = peakVal;
        end
    end
end

% ===== Sort by peak descending and select top N =====
[~, sortIdx] = sort([filePeaks.peak], 'descend');
topFiles = filePeaks(sortIdx(1:min(maxFiles, length(sortIdx))));

% ===== Plotting Loop =====
figure('Name','Aligned Spectra, Spectrograms & Waveforms','Color','w');

for fileCount = 1:length(topFiles)
    i = topFiles(fileCount).index;
    filePath = fullfile(files(i).folder, files(i).name);
    [data, fs] = audioread(filePath);

    % ===== ALIGN ALL CHANNELS TO CHANNEL 1 ======
    alignedData = data;
    ref = data(:,1);
    for ch = 2:4
        [xc, lags] = xcorr(data(:,ch), ref, 'coeff');
        [~, idx] = max(xc);
        lag = lags(idx);
        alignedData(:,ch) = circshift(data(:,ch), -lag);
    end

    % ===== Centre max peak of each channel ======
    N = size(alignedData, 1);
    midpoint = round(N / 2);
    for ch = 1:4
        [~, peakIdx] = max(abs(alignedData(:,ch)));
        shift = midpoint - peakIdx;
        alignedData(:,ch) = circshift(alignedData(:,ch), shift);
    end

    row = fileCount;

    % ====== Plot 1: pspec of aligned channels ======
    subplot(maxFiles, 4, (row - 1)*4 + 1);
    hold on;
    nfft = 1024;
    win = hann(nfft);
    noverlap = round(0.5 * nfft);
    f = (0:nfft/2) * fs / nfft;

    for ch = 1:4
        [Pxx, ~] = pwelch(alignedData(:,ch), win, noverlap, nfft, fs);
        Pxx = Pxx(1:nfft/2+1);
        Pxx = Pxx ./ max(Pxx);  % Normalise to [0,1]

        baseRGB = getColorRGB(chColors{ch});
        plot(f/1e3, 10*log10(Pxx), 'Color', baseRGB, 'LineWidth', 0.5);  % freq in kHz
    end
    xlim([0 fs/2000]); % display up to Nyquist in kHz
    ylim([-40 10]);
    xlim([20 80]);
    ylabel('dB, Norm.');

    if fileCount == 1
        title(sprintf('Normalised Spectrum'));
    end
    if fileCount == maxFiles
        xlabel('Frequency (kHz)');
    end
    formatLatex(gca)
    % grid on;

    % ====== Plot 2: img_spec on aligned channel 1 ======
    subplot(maxFiles, 4, (row - 1)*4 + 2);
    try
        img_spec(alignedData(:,1), 128, 127, 4096, fs, 100);
        if fileCount == 1
            title('Spectrogram');
        end
        yticklabels(yticks./1000)
        ylabel('kHz')
        if fileCount == maxFiles
            xlabel('Time (ms)');
        end
        formatLatex(gca)
    catch ME
        warning("img_spec failed on %s: %s", files(i).name, ME.message);
    end

    % ====== Plot 3: waveform of all aligned channels ======
    subplot(maxFiles, 4, (row - 1)*4 + 3);
    hold on;
    t = (0:length(alignedData)-1) ./ fs * 1000;  % time in ms
    for ch = 1:4
        % Convert base colour + alpha
        baseRGB = getColorRGB(chColors{ch});  % helper to get RGB triplet
        plot(t, alignedData(:,ch), 'Color', [baseRGB 0.2], 'LineWidth', 0.5);
    end
    if fileCount == 1
        title('Aligned Waveforms');
    end
    if fileCount == maxFiles
        xlabel('Time (ms)');
    end
    % ylabel('Amplitude');
    % xlim([min(t), max(t)]);
    xlim([5 13])
    formatLatex(gca)
    % grid on;


    % Prepare row entry for duration table
    durRow = struct('File', files(i).name, 'Dur1', NaN, 'Dur2', NaN, 'Dur3', NaN, 'Dur4', NaN);


    for ch = 1:4
        sig = alignedData(:,ch);

        % --- Use percentdur function with [5% 95%] cumulative energy ---
        [durationIdx, ~, ~] = percentdur(sig, [0.1 0.90]);

        % Convert indices to time
        leftTime = t(durationIdx(1));
        rightTime = t(durationIdx(2));

        % Plot vertical lines
        xline(leftTime, '--', 'Color', chColors{ch}, 'LineWidth', 0.5);
        xline(rightTime, '--', 'Color', chColors{ch}, 'LineWidth', 0.5);

        % Plot waveform
        plot(t, sig, 'Color', chColors{ch});

        % Save duration in ms
        durRow.(['Dur' num2str(ch)]) = rightTime - leftTime;
    end

    % title('Aligned Waveforms');
    % xlabel('Time (ms)');
    % ylabel('Amplitude');
    % xlim([min(t), max(t)]);
    xlim([5 13])
    % formatLatex(gca)


    % Add to duration table
    if fileCount == 1
        durationTable = struct2table(durRow);
    else
        durationTable = [durationTable; struct2table(durRow)];
    end

    % --- Compute AM envelopes for all 4 channels ---
    envs = zeros(size(alignedData));
    for ch = 1:4
        sig = alignedData(:,ch);
        envs(:,ch) = filtfilt(butter(4, 5000/(fs/2), 'low'), 1, abs(hilbert(sig)));
    end

    % ====== Plot 4: AM Envelopes (Normalised) ======
    subplot(maxFiles, 4, (row - 1)*4 + 4);
    hold on;

    for ch = 1:4
        baseRGB = getColorRGB(chColors{ch});
        envNorm = envs(:,ch) ./ max(envs(:,ch));  % Normalise to [0, 1]
        plot(t, envNorm, 'Color', [baseRGB 0.8], 'LineWidth', 0.5);
    end
    if fileCount == 1
        title('Normalised AM Envelopes');
    end
    if fileCount == maxFiles
        xlabel('Time (ms)');
    end
    ylabel('Norm. Envelope');
    xlim([5 13]);
    ylim([0 1.05]);
    formatLatex(gca);

    % --- Compare all unique pairs of channels ---
    pairIdx = 0;
    for ch1 = 1:3
        for ch2 = (ch1+1):4
            env1 = envs(:,ch1);
            env2 = envs(:,ch2);
            n = min(length(env1), length(env2));
            env1 = env1(1:n);
            env2 = env2(1:n);

            % Pearson Correlation
            r = corr(env1(:), env2(:));

            % RMSE
            rmse = sqrt(mean((env1 - env2).^2));

            % Cross-correlation
            [xc, lags] = xcorr(env1, env2, 'coeff');
            [max_corr, idx] = max(xc);
            lag_ms = lags(idx) / fs * 1000;

            % DTW (optional)
            try
                dist_dtw = dtw(env1, env2);
            catch
                dist_dtw = NaN;
            end

            % Store each pair's result
            pairIdx = pairIdx + 1;
            compRow = struct( ...
                'File', fileCount, ...
                'Channel1', ch1, ...
                'Channel2', ch2, ...
                'PearsonCorr', str2double(sprintf('%.2e', r)), ...
                'RMSE',        str2double(sprintf('%.2e', rmse)), ...
                'MaxXCorr',    str2double(sprintf('%.2e', max_corr)), ...
                'LagMS',       str2double(sprintf('%.2e', lag_ms)), ...
                'DTW',         str2double(sprintf('%.2e', dist_dtw)));

            if fileCount == 1 && pairIdx == 1
                amSummaryTable = struct2table(compRow);
            else
                amSummaryTable = [amSummaryTable; struct2table(compRow)];
            end
        end
    end
end
sgtitle('Cross-Channel Signal Characterisation - Spatial Invariance', 'FontSize', 18, 'FontWeight', 'bold', 'Interpreter', 'latex')
%% save
% saveFigure(gcf, fig_path, 'spatial_invariance')

% writetable(amSummaryTable, fullfile('../results', 'am_summary_table.csv'));

%%
function rgb = getColorRGB(colorChar)
switch colorChar
    case 'r', rgb = [1, 0, 0];
    case 'g', rgb = [0, 1, 0];
    case 'b', rgb = [0, 0, 1];
    case 'k', rgb = [0, 0, 0];
    case 'm', rgb = [1, 0, 1];
    case 'c', rgb = [0, 1, 1];
    case 'y', rgb = [1, 1, 0];
    case 'w', rgb = [1, 1, 1];
    otherwise, rgb = [0.5, 0.5, 0.5];  % default grey
end
end