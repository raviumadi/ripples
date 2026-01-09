classdef RippleStudio < handle
    % RIPPLESTUDIO V1.2  Interactive demo of spectral interference ("ripples")
    % in bat echolocation over water surfaces.
    %
    % RippleStudio simulates a down-swept FM bat call emitted above a water
    % surface. The direct path and its water-reflected copy interfere,
    % producing the characteristic spectral ripple (an AM signature in the
    % frequency domain). The app lets you explore how geometry and call
    % parameters shape this interference, with three linked visualisations:
    %   • Spectrogram of the overlapped (direct + delayed) call
    %   • Time waveform of the overlapped call
    %   • Plan-view geometry with optional beam boundary lines
    %
    % =========================================================================
    % QUICK START
    %   h = RippleStudio;   % launch GUI
    %
    %   % Export any axes exactly “as seen” (titles, labels, ticks preserved):
    %   h.exportAxestoPDF(h.ax1,  'spectrogram.pdf');   % spectrogram
    %   h.exportAxestoPDF(h.ax1b, 'waveform.pdf');      % waveform
    %   h.exportAxestoPDF(h.ax2,  'diagram.pdf');       % geometry
    %
    %   % Programmatic access to generated signals:
    %   x      = h.call;            % direct FM call (column)
    %   x_ref  = h.indirectCall;    % delayed/reflected copy
    %   x_sum  = h.overlapCall;     % 0.5*(x + x_ref)
    %   t_ms   = h.t_ms;            % time vector (ms)
    %   dt_s   = h.dt_s;            % geometric delay τ (s)
    %
    % =========================================================================
    % CONTROLS
    %   Sliders (left column)
    %     • Bat Height (m)            – emitter altitude above water
    %     • Beam Angle (°)            – downward angle of sonar axis
    %     • Shore Distance (m)        – horizontal offset to shoreline
    %     • Call Duration (s)         – FM sweep length
    %     • Aperture Diameter (m)     – piston source size (beamwidth)
    %     • Start Freq (kHz)          – FM sweep start (high) frequency
    %     • End Freq (kHz)            – FM sweep end (low) frequency
    %     • Beam Boundary (dB)        – cutoff for beam edges (e.g., −3 dB)
    %
    %   Checkboxes
    %     • Show start-freq. boundary – beam edges at Start Freq
    %     • Show end-freq. boundary   – beam edges at End Freq
    %
    %   Buttons
    %     • Reset     – restore all sliders to defaults
    %     • ⓘ / Help  – open scrollable documentation window
    %
    % =========================================================================
    % VISUAL OUTPUTS
    %   Spectrogram (top-left):  Overlapped call. Title reports θ, h, and f_rpl,
    %     where τ = (L_reflected − L_direct)/c, f_rpl = 1/τ (Hz).
    %     Null-to-null spectral spacing equals f_rpl; temporal beat is
    %     f_AM = 1/(2τ).
    %
    %   Waveform (top-right):    Overlapped signal in time.
    %
    %   Geometry (bottom):       Direct/reflected paths, interference point,
    %                            and optional beam edges at selected cutoff.
    %
    % =========================================================================
    % SIGNALS (public get)
    %   h.call           – direct FM call
    %   h.indirectCall   – delayed/reflected copy
    %   h.overlapCall    – 0.5*(call + indirectCall)
    %   h.t_ms           – time vector (ms)
    %   h.dt_s           – geometric delay τ (s)
    %
    % =========================================================================
    % METHODS (public)
    %   exportAxestoPDF(ax, filename, ...)
    %       Export any axes “as seen” (preserves formatting).
    %       Defaults to PDF, image content, 300 dpi.
    %       Examples:
    %         h.exportAxestoPDF(h.ax1, 'spec.pdf');
    %         h.exportAxestoPDF(h.ax2, 'diagram.pdf', 'Resolution', 600);
    %
    % =========================================================================
    % MODEL
    %   • Call: Down-swept FM (quadratic), Hann-tapered, optional padding.
    %   • Beamwidth: Circular piston model; edges at chosen dB cutoff.
    %   • Interference: Overlap of direct and delayed copies (mean).
    %   • Spectrogram: |STFT| in dB, fixed dynamic range.
    %
    % =========================================================================
    % MATH REFERENCE
    %   c      = 343 m/s
    %   τ      = (L_reflected − L_direct)/c   (s)   geometric delay
    %   f_rpl  = 1/τ                          (Hz)  spectral ripple spacing
    %   f_AM   = 1/(2τ)                       (Hz)  temporal beat
    %
    % =========================================================================
    % FILES & REQUIREMENTS
    %   • documentation.txt (optional): populates Help panel.
    %   • MATLAB R2018b+ (R2021a+ recommended for exportgraphics).
    %
    % =========================================================================
    % CITATION
    %   If used in publications or teaching, please cite:
    %   Umadi & Firzlaff (2025) and acknowledge RippleStudio.
    %
    % =========================================================================
    % AUTHOR / LICENSE
    %   Author: <>
    %   License: CC-BY-NC 4.0. See LICENSE file.
    % -------------------------------------------------------------------------
    properties (Constant)
        c  = 343;          % Speed of sound (m/s)
        fs = 192e3;        % Sampling rate
    end

    properties
        f
        ax1, ax1b, ax2
        infoBox
        infoText
        cbStart
        cbEnd
        infoAnn   % LaTeX-capable bottom info box
        helpBtn
        helpFig
        helpAx
        helpSlider
        helpTextHandles
        resetBtn
        sliderHandles    % handles to all sliders (creation order)
        sliderDefaults   % default values parallel to sliderHandles
    end

    properties (GetAccess = public, SetAccess = private)
        call            % direct FM call (column)
        indirectCall    % delayed/reflected copy (column)
        overlapCall     % 0.5 * (call + indirectCall)
        t_ms            % time vector in ms (column)
        dt_s            % delay between paths (seconds)
        param           % all calculated parameters
    end

    methods
        function self = RippleStudio()
            % --- GUI Figure (keep layout style) ---
            self.f = figure('Name', "Ripple Studio V1.2 - Ravi Umadi (2025)", ...
                'NumberTitle', 'off', 'Position', [100 100 1250 750], 'Color', [0.96 1 0.96], 'Resize', 'off'); % [0.96 1 0.96]

            % --- Axes (same positions as your original) ---
            self.ax1  = axes('Parent', self.f, 'Position', [0.35 0.55 0.28 0.4]); % spectrogram
            self.ax1b = axes('Parent', self.f, 'Position', [0.67 0.55 0.28 0.4]); % waveform
            self.ax2  = axes('Parent', self.f, 'Position', [0.35 0.15 0.6 0.3]);  % diagram
            set([self.ax1,self.ax1b,self.ax2], 'TickLabelInterpreter','latex', 'FontSize', 12, 'Box','on');

            % --- UI Controls (labels with LaTeX) ---
            % ----- Sliders (auto-stacked) -----
            specs = {
                '$\mathbf{Bat\ Height}\ (\mathrm{m})$',              0.1,    1.5,   0.2,   []
                '$\mathbf{Beam\ Angle}\ (^\circ)$',                  1,      89,    28,    []
                '$\mathbf{Shore\ Distance}\ (\mathrm{m})$',          0.5,    25,    8,     []
                '$\mathbf{Call\ Duration}\ (\mathrm{s})$',           0.002,  0.01,  0.004, []
                '$\mathbf{Aperture\ Diameter}\ (\mathrm{m})$',       0.005,  0.02,  0.018, []
                '$\mathbf{Start\ Freq}\ (\mathrm{kHz})$',            60,     95,    90,    []
                '$\mathbf{End\ Freq}\ (\mathrm{kHz})$',              25,     50,    25,    []
                '$\mathbf{Beam\ Boundary}\ (\mathrm{dB})$',          -12,    -3,    -3,    []
                '$\mathbf{Responsivity\ Coeff.}\ (k_r)$',            1,      10,    5,     []
                };

            baseY  = 660;   % top slider y-position (pixels)
            stepY  = 60;    % vertical spacing (≥85 recommended for label + slider + ticks)
            leftX  = 50;    % left x (pixels)
            width  = 250;   % slider width
            height = 15;    % slider height

            self.addSliderStack(specs, baseY, stepY, leftX, width, height, @self.updatePlot);

            % --- Reset button (bottom of slider stack) ---
            nSliders = size(specs,1);
            bottomY  = baseY - (nSliders-1)*stepY;     % y of the lowest slider
            btnY     = bottomY - 110;                   % a bit below the lowest slider
            btnH     = 24;
            self.resetBtn = uicontrol(self.f, 'Style','pushbutton', ...
                'String','Reset', ...
                'Position',[leftX btnY width btnH], ...
                'TooltipString','Restore all sliders to their defaults', ...
                'Callback', @(~,~) self.resetSliders());

            % === Checkboxes with LaTeX labels ===
            % Start boundary
            self.cbStart = uicontrol(self.f, 'Style','checkbox', ...
                'Value',1, 'Position',[50 130 20 20], 'Callback', @self.updatePlot);

            figPos = getpixelposition(self.f);   % figure size in pixels
            labPix = [75 130 250 20];            % label rectangle in pixels (just right of box)
            labNorm = [labPix(1:2)./figPos(3:4), labPix(3:4)./figPos(3:4)];

            annotation(self.f,'textbox',labNorm, ...
                'String','$\mathrm{Show\ start\!-\!freq.\ boundary}$', ...
                'Interpreter','latex', ...
                'EdgeColor','none','BackgroundColor','none', ...
                'HorizontalAlignment','left','VerticalAlignment','middle', ...
                'FontSize',12,'Margin',0);

            % End boundary
            self.cbEnd = uicontrol(self.f, 'Style','checkbox', ...
                'Value',1, 'Position',[50 100 20 20], 'Callback', @self.updatePlot);

            labPix2 = [75 100 250 20];
            labNorm2 = [labPix2(1:2)./figPos(3:4), labPix2(3:4)./figPos(3:4)];

            annotation(self.f,'textbox',labNorm2, ...
                'String','$\mathrm{Show\ end\!-\!freq.\ boundary}$', ...
                'Interpreter','latex', ...
                'EdgeColor','none','BackgroundColor','none', ...
                'HorizontalAlignment','left','VerticalAlignment','middle', ...
                'FontSize',12,'Margin',0);

            % --- Info Displays (bottom) ---
            % --- Bottom info: single LaTeX annotation box (pixel units)
            posPix = [50 20 1150 40];  % [x y w h] in pixels (w = fixed width)
            self.infoAnn = annotation(self.f,'textbox',[0 0 1 1], ...
                'Interpreter','latex', 'String','', ...
                'FontSize',11, ...
                'EdgeColor',[0.8 0.8 0.8], 'BackgroundColor',[1 0.95 0.95], ...
                'HorizontalAlignment','left','VerticalAlignment','middle', ...
                'Margin',8,'FitBoxToText','on','LineWidth',0.5);
            set(self.infoAnn,'Units','pixels','Position',posPix);

            % ----- Help menu -----
            mHelp = uimenu(self.f,'Label','Help');
            uimenu(mHelp,'Label','Show description…','Accelerator','H', ...
                'Callback', @(~,~) self.showDescription());

            % ----- Tiny "info" button (top-right) -----
            self.helpBtn = uicontrol(self.f,'Style','pushbutton', ...
                'String','ⓘ','FontWeight','bold','FontSize',12, ...
                'TooltipString','Description', ...
                'BackgroundColor',[0.94 0.97 1.00], ...
                'Callback', @(~,~) self.showDescription());

            % place it initially & keep it in the corner on resize
            self.placeHelpButton();                      % initial placement
            self.f.SizeChangedFcn = @(~,~) self.placeHelpButton();  % keep pinned

            % First draw
            self.updatePlot();
        end
    end

    methods (Access = private)
        function addSlider(self, label, minVal, maxVal, initVal, pos, callback, steps)
            % --- LaTeX label above slider ---
            figPos  = getpixelposition(self.f);                 % [x y w h] pixels
            labPix  = [pos(1) pos(2)+15 260 20];                % label rectangle
            labNorm = [labPix(1:2)./figPos(3:4), labPix(3:4)./figPos(3:4)];
            annotation(self.f,'textbox', labNorm, ...
                'String', label, 'Interpreter','latex', ...
                'EdgeColor','none','BackgroundColor','none', ...
                'HorizontalAlignment','left','VerticalAlignment','middle', ...
                'FontSize', 12, 'Margin',0);

            % --- Slider (slightly styled) ---
            s = uicontrol(self.f,'Style','slider', ...
                'Min',minVal,'Max',maxVal,'Value',initVal, ...
                'Position',pos,'Callback',callback, ...
                'BackgroundColor',[0.92 0.93 0.96]); % subtle tint
            if nargin == 8 && ~isempty(steps), set(s,'SliderStep',steps); end
            s.Tag = label;

            % keep a handle + its default for resetting later
            if isempty(self.sliderHandles), self.sliderHandles = gobjects(0); end
            self.sliderHandles(end+1)  = s;
            self.sliderDefaults(end+1) = initVal;

            % --- Fancy ticks BELOW the slider ---
            nticks   = 5;
            tickVals = linspace(minVal,maxVal,nticks);

            % choose decimals based on step size (prevents 0.01,0.01,...)
            step = (maxVal - minVal) / max(1, nticks-1);
            if step < 1e-3
                fmt = '%.5f';
            elseif step < 1e-2
                fmt = '%.4f';
            elseif step < 1e-1
                fmt = '%.3f';
            else
                fmt = '%.2f';
            end

            tickStr  = arrayfun(@(v) num2str(v, fmt), tickVals, 'UniformOutput', false);

            tickHeight = 12;
            tickOffset = 10;
            tickPix    = [pos(1) pos(2)-tickOffset pos(3) tickHeight];
            tickNorm   = [tickPix(1:2)./figPos(3:4), tickPix(3:4)./figPos(3:4)];
            axTicks = axes('Parent', self.f, 'Position', tickNorm, 'Color','none', ...
                'XLim',[minVal maxVal], 'YLim',[0 1], 'XTick',tickVals, ...
                'XTickLabel',tickStr, 'YTick',[], 'Box','off', ...
                'TickDir','out', 'FontSize', 9);
            axTicks.XAxis.TickLabelInterpreter = 'latex';
            uistack(axTicks,'bottom');

            % --- Right-end numeric readout (use same fmt) ---
            readoutPix  = [pos(1)+pos(3)+8, pos(2)-2, 64, 20];
            readoutNorm = [readoutPix(1:2)./figPos(3:4), readoutPix(3:4)./figPos(3:4)];
            valTxt = annotation(self.f,'textbox', readoutNorm, ...
                'String', num2str(initVal, fmt), ...
                'Interpreter','latex','EdgeColor','none', ...
                'FontSize', 12, 'HorizontalAlignment','left', 'VerticalAlignment','middle');

            addlistener(s,'Value','PostSet', @(~,evt) ...
                set(valTxt,'String', num2str(evt.AffectedObject.Value, fmt)));

        end

        function positionValueLabel(~, txt, slider, pos, minVal, maxVal, val)
            % compute relative x offset
            frac = (val-minVal)/(maxVal-minVal);
            x = pos(1) + frac*pos(3);
            y = pos(2) + pos(4) + 5;
            set(txt,'Units','pixels','Position',[x-15 y 30 20], ...
                'String',num2str(val,'%.2f'));
        end

        function vals = getSliders(self)
            % Fetch all sliders sorted by vertical position (top -> bottom)
            all = findall(self.f, 'Style', 'slider');
            [~, idx] = sort(arrayfun(@(h) h.Position(2), all), 'descend');
            all = all(idx);
            vals = arrayfun(@(h) h.Value, all);
        end

        function resetSliders(self)
            % Restore each slider to its original default value
            if isempty(self.sliderHandles) || isempty(self.sliderDefaults), return; end
            for k = 1:numel(self.sliderHandles)
                h = self.sliderHandles(k);
                if isgraphics(h) && isvalid(h)
                    set(h, 'Value', self.sliderDefaults(k));   % triggers PostSet -> readout updates
                end
            end
            % one redraw at the end
            self.updatePlot();
        end

        function addSliderStack(self, specs, baseY, stepY, leftX, width, height, callback)
            % specs: N×5 cell array:
            %   {label, min, max, init, stepsOrEmpty}
            % Spacing:
            %   baseY = top y (pixels), stepY = vertical spacing (pixels)
            %   leftX/width/height = slider rect (pixels)
            %
            % Tip: with the “fancy ticks” addSlider, use stepY >= 85 (label+slider+ticks).

            for i = 1:size(specs,1)
                label   = specs{i,1};
                minVal  = specs{i,2};
                maxVal  = specs{i,3};
                initVal = specs{i,4};
                if size(specs,2) >= 5 && ~isempty(specs{i,5})
                    steps = specs{i,5};
                else
                    steps = [];
                end
                y = baseY - (i-1)*stepY;
                self.addSlider(label, minVal, maxVal, initVal, [leftX y width height], callback, steps);
            end
        end

        function updatePlot(self, ~, ~)
            % === Gather control values (order = by Y position) ===
            vals = self.getSliders();
            % Expected order (top->bottom): h, theta, shore, duration, aperture, f_start(kHz), f_end(kHz), cutoff(dB)
            h            = vals(1);
            theta_deg    = vals(2);
            shore_dist   = vals(3);
            call_duration= vals(4);
            aperture     = vals(5);
            f_start_kHz  = vals(6);
            f_end_kHz    = vals(7);
            cutoff_db    = vals(8);
            kr          = vals(9);   % responsivity coefficient from slider

            f_start = f_start_kHz*1e3;
            f_end   = f_end_kHz*1e3;

            % === Generate call ===
            self.call = self.generateVirtualBatCall(f_start, f_end, call_duration, self.fs, 100);
            theta = deg2rad(theta_deg);

            % === Geometry (your original) ===
            bat = [0, h];
            x_r = h / tan(theta);
            reflection   = [x_r, 0];
            dx = h / tan(theta);
            interference = reflection + [dx, h];
            shore = [shore_dist, h];

            a = norm(bat - reflection);
            b = norm(reflection - interference);
            direct = norm(bat - interference);
            total_reflected = a + b;
            path_diff = total_reflected - direct;
            delay_ms = 1000 * path_diff / self.c;
            tau = delay_ms;
            f_ripple = 1/(tau/1000);

            % Beam widths at cutoff_db
            beam_start = self.computePistonSonarBeamWidth(aperture, f_start, cutoff_db);
            beam_end   = self.computePistonSonarBeamWidth(aperture, f_end,   cutoff_db);

            % === Call rate from beamwidth-derived distance ===
            beam_len = 2*h * tan(deg2rad(beam_end.beamwidth_deg/2));  % lateral footprint
            d_eff    = sqrt(h^2 + (beam_len/2)^2);                    % effective slant range
            Cr_eff   = self.c ./ (2 * (1 + kr) * d_eff);              % responsivity model

            % === Two-path interference in time ===
            dt = path_diff / self.c;
            indirect_call = circshift(self.call, round(self.fs * dt));
            overlap_call = mean([self.call indirect_call], 2);

            % === Save to object for external access ===
            self.indirectCall = indirect_call(:);
            self.overlapCall  = overlap_call(:);
            self.dt_s         = dt;
            self.t_ms         = (0:numel(overlap_call)-1)'/self.fs * 1000;

            % === Spectrogram (img_spec) ===
            axes(self.ax1); cla(self.ax1);
            self.img_spec(overlap_call, 128, 127, 4096, self.fs, 100);
            grid(self.ax1,'on');
            ylim(self.ax1,[0 100e3]);
            yt = yticks(self.ax1); yticklabels(self.ax1, yt./1000); % kHz
            xt = xticks(self.ax1); xticklabels(self.ax1, xt.*1000); % ms
            title(self.ax1, ...
                sprintf('\\textbf{Spectrogram} ($\\theta = %.1f^{\\circ},\\ h = %.2f\\,\\mathrm{m},\\ f_{\\mathrm{rpl}} = %.2f\\,\\mathrm{kHz}$)', ...
                theta_deg, h, f_ripple/1000), ...
                'Interpreter','latex');
            xlabel(self.ax1,'$t\ (\mathrm{ms})$','Interpreter','latex');
            ylabel(self.ax1,'$f\ (\mathrm{kHz})$','Interpreter','latex');
            set(self.ax1, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

            % === Waveform ===
            axes(self.ax1b); cla(self.ax1b);
            t = (0:length(overlap_call)-1)/self.fs * 1000;
            plot(self.ax1b, t, overlap_call, 'k');
            xlabel(self.ax1b,'$t\ (\mathrm{ms})$','Interpreter','latex');
            ylabel(self.ax1b,'$\mathrm{Amplitude}$','Interpreter','latex');
            title(self.ax1b, sprintf('\\textbf{Overlapped Call}'), ...
                'Interpreter','latex');
            grid(self.ax1b,'on'); xlim(self.ax1b,[0 max(t)]);
            set(self.ax1b, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

            % === Diagram (with boundary toggles) ===
            axes(self.ax2); cla(self.ax2); hold(self.ax2,'on');
            plot(self.ax2, [bat(1), reflection(1)], [bat(2), reflection(2)], 'b-', 'LineWidth', 2);
            plot(self.ax2, [reflection(1), interference(1)], [reflection(2), interference(2)], 'r-', 'LineWidth', 2);
            plot(self.ax2, [bat(1), interference(1)], [bat(2), interference(2)], 'g--', 'LineWidth', 1.5);
            plot(self.ax2, [bat(1), shore(1)], [h h], 'k--', 'LineWidth', 1);
            yl = yline(0, 'k', 'Water', 'LineWidth',1);
            yl.LabelHorizontalAlignment = 'center';   % 'left' | 'center' | 'right'
            yl.LabelVerticalAlignment   = 'bottom';    % 'top'  | 'middle' | 'bottom'
            yl.FontSize = 12;
            yl.Interpreter = 'latex';               % if you want LaTeX formatting
            plot(self.ax2, bat(1), bat(2), 'ko', 'MarkerFaceColor', 'b');
            plot(self.ax2, reflection(1), reflection(2), 'ko', 'MarkerFaceColor', 'r');
            plot(self.ax2, interference(1), interference(2), 'ko', 'MarkerFaceColor', 'm');
            plot(self.ax2, shore(1), shore(2), 'ko', 'MarkerFaceColor', 'k');

            text(self.ax2, shore(1), shore(2)+0.1, 'Shore', 'HorizontalAlignment', 'center', 'Fontsize', 12, 'FontWeight', 'bold', 'interpreter', 'latex');
            text(self.ax2, bat(1), bat(2)+0.1, 'Bat', 'HorizontalAlignment', 'center', 'Fontsize', 12, 'FontWeight', 'bold', 'interpreter', 'latex');
            text(self.ax2, interference(1), interference(2)+0.1, 'I', 'HorizontalAlignment', 'center', 'Fontsize', 12, 'FontWeight', 'bold', 'interpreter', 'latex');
            text(self.ax2, (bat(1)+reflection(1))/2, h/2,  'a' , 'Color', 'b', 'HorizontalAlignment', 'left', 'Fontsize', 12, 'FontWeight', 'bold', 'interpreter', 'latex');
            text(self.ax2, (reflection(1)+interference(1))/2, h/2, 'b', 'Color', 'r', 'HorizontalAlignment', 'left', 'Fontsize', 12, 'FontWeight', 'bold', 'interpreter', 'latex');

            % Show start/end beam boundaries only if enabled
            if self.cbStart.Value
                self.drawBeamEdges(bat, h, theta, beam_start, [1, 0.5, 0]); % orange
            end
            if self.cbEnd.Value
                self.drawBeamEdges(bat, h, theta, beam_end,   [0, 0.5, 0.5]); % teal
            end

            xlim(self.ax2, [-0.5, max(shore_dist+1, interference(1)+1)]);
            ylim(self.ax2, [-0.2, h+1]);
            xlabel(self.ax2,'$x\ (\mathrm{m})$','Interpreter','latex');
            ylabel(self.ax2,'$y\ (\mathrm{m})$','Interpreter','latex');
            % First line (text mode)
            t1 = sprintf('\\textbf{Reflection Diagram} ($\\theta = %.1f^{\\circ}$)', theta_deg);

            % Second line (math mode) — compact, single string
            longtitle = {
                ['\mathbf{a,b}=', num2str(a,'%.2f'), '\,\mathrm{m}']
                ['\mathbf{C_{r}}=', num2str(Cr_eff,'%.2f'), '\,\mathrm{Hz}']
                % ['\mathbf{b}=', num2str(b,'%.2f'), '\,\mathrm{m}']
                % ['\mathbf{L_{\mathrm{refl}}}=', num2str(total_reflected,'%.2f'), '\,\mathrm{m}']
                % ['\mathbf{L_{\mathrm{dir}}}=',  num2str(direct,'%.2f'),         '\,\mathrm{m}']
                ['\mathbf{\Delta L}=',          num2str(path_diff,'%.2f'),      '\,\mathrm{m}']
                ['\mathbf{\Delta t}=',          num2str(delay_ms,'%.2f'),       '\,\mathrm{ms}']
                ['\mathbf{D}=',                 num2str(aperture*1000,'%.1f'),  '\,\mathrm{mm}']
                ['\mathbf{BW}=', num2str(beam_start.beamwidth_deg,'%.2f'), '^{\circ}', '/' num2str(beam_end.beamwidth_deg,'%.2f'), '^{\circ}']
                % ['\mathbf{BW_{end}}=', num2str(beam_end.beamwidth_deg,'%.2f'),   '^{\circ}']
                ['\mathbf{Cutoff}=',    num2str(cutoff_db,'%.0f'),     '\,\mathrm{dB}']
                };
            t2 = ['$', strjoin(longtitle, '\\quad '), '$'];  % single LaTeX line

            % Two-line title: main title + compact subtitle
            title(self.ax2, {t1, t2}, 'Interpreter','latex', 'FontSize', 11);
            grid(self.ax2,'on');

            % === Info (bottom) — single line, LaTeX (no \text{...}) ===
            parts = {
                ['a, b=',           num2str(a,'%.2f'),           '\,\mathrm{m}']
                ['f_{rpl}=',           num2str(f_ripple/1000,'%.2f'),           '\,\mathrm{kHz}']
                ['L_{\mathrm{refl}}=', num2str(total_reflected,'%.2f'), '\,\mathrm{m}']
                ['L_{\mathrm{dir}}=',  num2str(direct,'%.2f'),         '\,\mathrm{m}']
                ['\Delta L=',    num2str(path_diff,'%.2f'),   '\,\mathrm{m}']
                ['\Delta t=',    num2str(delay_ms,'%.2f'),    '\,\mathrm{ms}']
                ['h=',           num2str(h,'%.2f'),           '\,\mathrm{m}']
                ['\theta=',      num2str(theta_deg,'%.1f'),   '^{\circ}']
                ['\mathrm{Shore}=', num2str(shore_dist,'%.2f'), '\,\mathrm{m}']
                ['T=',           num2str((call_duration*1000),'%.1f'), '\,\mathrm{ms}']
                ['D=',           num2str((aperture*1000),'%.1f'),      '\,\mathrm{mm}']
                ['f_{\mathrm{start}}=', num2str(f_start/1e3,'%.1f'),  '\,\mathrm{kHz}\Rightarrow \mathrm{BW}=', num2str(beam_start.beamwidth_deg,'%.2f'), '^{\circ}']
                ['f_{\mathrm{end}}=',   num2str(f_end/1e3,'%.1f'),    '\,\mathrm{kHz}\Rightarrow \mathrm{BW}=', num2str(beam_end.beamwidth_deg,'%.2f'),   '^{\circ}']
                ['\mathrm{Cut}=',    num2str(cutoff_db,'%.0f'),    '\,\mathrm{dB}']
                };

            % --- Choose how many chunks per line (tune for your FontSize/width)
            perLine = 2;  % try 5–6 for FontSize ~18 and width ~1100 px

            % --- Build final LaTeX string and set it
            infoStr = self.buildInfoLatex(parts, perLine);
            set(self.infoAnn,'String', infoStr);

            % === Collect all parameters into a struct ===
            self.param = struct( ...
                'h', h, ...
                'theta_deg', theta_deg, ...
                'shore_dist', shore_dist, ...
                'call_duration_ms', call_duration*1000, ...
                'aperture_mm', aperture*1000, ...
                'f_start_kHz', f_start/1e3, ...
                'f_end_kHz', f_end/1e3, ...
                'cutoff_db', cutoff_db, ...
                'kr', kr, ...
                'a', a, ...
                'b', b, ...
                'L_refl', total_reflected, ...
                'L_dir', direct, ...
                'deltaL', path_diff, ...
                'deltaT_ms', delay_ms, ...
                'f_rpl_kHz', f_ripple/1000, ...
                'Cr_Hz', Cr_eff, ...
                'BW_start_deg', beam_start.beamwidth_deg, ...
                'BW_end_deg', beam_end.beamwidth_deg );

        end

        function drawBeamEdges(self, bat, h, theta, beam, color)
            if isnan(beam.beamwidth_deg), return; end
            half_angle = deg2rad(beam.beamwidth_deg / 2);
            len = 1.0;
            for sign = [-1, 1]
                phi = theta + sign * half_angle;
                line_width = (sign < 0) * 1 + (sign > 0) * 0.75;
                if phi < 0
                    x1 = bat(1) + len * cos(2*pi - phi);
                    y1 = bat(2) + len * sin(2*pi - phi);
                    plot(self.ax2, [bat(1), x1], [bat(2), y1], ':', 'Color', color, 'LineWidth', line_width);
                else
                    x_beam = h / tan(phi);
                    plot(self.ax2, [bat(1), x_beam], [bat(2), 0], '-', 'Color', color, 'LineWidth', line_width);
                    x_reflected = 2 * x_beam;
                    plot(self.ax2, [x_beam, x_reflected], [0, h], '--', 'Color', color, 'LineWidth', line_width);
                end
            end
        end

        function beam = computePistonSonarBeamWidth(~, apertureDiameter, frequency, cutoff_db)
            % (Your original helper, unmodified except cutoff param)
            c = RippleStudio.c;
            lambda = c / frequency;
            k = 2 * pi / lambda;
            a = apertureDiameter / 2;

            theta = linspace(-pi/2, pi/2, 5000);
            x = k * a * sin(theta);
            x(x == 0) = eps;

            pattern = abs(2 * besselj(1, x) ./ x);
            pattern = pattern / max(pattern);
            pattern(pattern < eps) = eps;

            pattern_dB = 20 * log10(pattern);

            [~, peakIdx] = max(pattern_dB);
            leftIdx  = find(pattern_dB(1:peakIdx) < cutoff_db, 1, 'last');
            rightIdx = find(pattern_dB(peakIdx:end) < cutoff_db, 1, 'first') + peakIdx - 1;

            if ~isempty(leftIdx) && ~isempty(rightIdx)
                theta_left = interp1(pattern_dB(leftIdx:leftIdx+1), theta(leftIdx:leftIdx+1), cutoff_db);
                theta_right = interp1(pattern_dB(rightIdx-1:rightIdx), theta(rightIdx-1:rightIdx), cutoff_db);
                beamwidth_deg = rad2deg(theta_right - theta_left);
            else
                beamwidth_deg = NaN;
            end

            beam = struct('frequency', frequency, 'apertureDiameter', apertureDiameter, ...
                'theta', theta, 'pattern_dB', pattern_dB, 'beamwidth_deg', beamwidth_deg);
        end

        function img_spec(~, sig, n, o, F, fs, range)
            % IMG_SPEC  Spectrogram image with fixed dynamic range.
            % Originally written by Lasse Jakobsen
            % Usage: self.img_spec(sig, n, o, F, fs, range)
            %
            % sig   : signal (column or row)
            % n     : window size
            % o     : overlap
            % F     : FFT length
            % fs    : sample rate (Hz)
            % range : linear dynamic range divisor (e.g., 100)

            % ensure column vector
            if size(sig,1) == 1
                sig = sig(:);
            end

            [B,f,t] = spectrogram(sig, n, o, F, fs);
            bmin    = max(max(abs(B))) / range;
            imagesc(t, f, 20*log10(max(abs(B), bmin) / bmin));
            ax = gca;
            set(ax, 'YDir','normal');
            colormap(ax, flipud(hot));
        end

        function s = buildInfoLatex(~, items, perLine)
            % items: cell array of small LaTeX chunks (no trailing \quad)
            % perLine: how many chunks per line before inserting a line break
            lines = {};
            for k = 1:perLine:numel(items)
                r = k:min(k+perLine-1, numel(items));
                lines{end+1} = strjoin(items(r), '\\quad '); %#ok<AGROW>
            end
            s = ['$', strjoin(lines, ' \\ '), '$'];  % wrap in math mode
        end

        function showDescription(self)
            % --- read doc file ---
            classPath = fileparts(mfilename('fullpath'));
            docFile   = fullfile(classPath, 'documentation.txt');
            if ~isfile(docFile)
                lines = {'$\mathbf{documentation.txt\ not\ found.}$'};
            else
                raw   = fileread(docFile);
                lines = strsplit(raw, newline);
                lines = lines(~cellfun(@isempty, strtrim(lines))); % drop empty edges
            end

            % --- (re)create window ---
            if ~isempty(self.helpFig) && isvalid(self.helpFig), close(self.helpFig); end
            self.helpFig = figure('Name','Documentation','NumberTitle','off', ...
                'Color',[1 0.95 0.95],'MenuBar','none','ToolBar','none', ...
                'WindowStyle','modal','Resize','off','Position',[100 100 700 800]);

            % Layout constants
            padL = 60; padR = 40; padT = 40; padB = 50;   % pixels
            sliderW = 16;
            fontSz = 16;
            lineStep = 0.075;   % vertical spacing (axes data units page=1)

            % --- axes for text (data-space scroll) ---
            axPos = self.localAxPos(self.helpFig, padL, padR+sliderW, padT, padB);
            self.helpAx = axes('Parent', self.helpFig, 'Position', axPos, ...
                'XLim',[0 1], 'YLim',[0 1], 'Visible','off');

            % --- draw lines once in data coords; we’ll scroll YLim later ---
            n = numel(lines);
            contentH = max(1, 0.04 + n*lineStep);   % content height in "pages"
            y = contentH - 0.03;                    % start near top
            self.helpTextHandles = gobjects(1,n);
            for k = 1:n
                L = strtrim(lines{k});
                if isempty(L), y = y - lineStep; continue; end
                self.helpTextHandles(k) = text(0.02, y, L, 'Parent', self.helpAx, ...
                    'Interpreter','latex', 'FontSize', fontSz, ...
                    'VerticalAlignment','top', 'HorizontalAlignment','left', ...
                    'Clipping','on');
                y = y - lineStep;
            end

            % --- vertical slider to scroll 1 page ---
            maxScroll = max(0, contentH - 1);  % page height = 1
            val0 = maxScroll;                  % show top initially

            % decide steps safely
            if maxScroll <= 0
                smallStep = 0;    % no scrolling needed
                largeStep = 0;
                sliderEnable = 'off';
            else
                % one "tick" ~ 1/numberOfLines per page, but cap nicely
                smallStep = min(0.25, 1/max(4, contentH));  % ≲ 0.25 page
                largeStep = max(smallStep, 0.2);            % ≥ smallStep
                sliderEnable = 'on';
            end

            self.helpSlider = uicontrol(self.helpFig,'Style','slider', ...
                'Units','pixels', ...
                'Position', self.localSliderPos(self.helpFig, sliderW, padR, padT, padB), ...
                'Min',0,'Max',maxScroll,'Value',val0, ...
                'Enable', sliderEnable, ...
                'SliderStep',[smallStep, largeStep], ...
                'Callback', @(s,~) set(self.helpAx,'YLim',[s.Value s.Value+1]));

            % initialize view
            set(self.helpAx,'YLim',[val0 val0+1]);

            % --- Close button pinned bottom-right ---
            btn = uicontrol(self.helpFig,'Style','pushbutton','String','Close', ...
                'Position',[self.helpFig.Position(3)-90, 10, 80, 28], ...
                'Callback', @(~,~) close(self.helpFig));

            % --- handle figure resize: relayout axes, slider, button ---
            self.helpFig.SizeChangedFcn = @(~,~) self.localOnHelpResize( ...
                padL,padR,padT,padB,sliderW,btn,contentH);
        end
        function pos = localAxPos(self, fig, padL, padR, padT, padB)
            old = fig.Units; fig.Units = 'pixels';
            fp = fig.Position; fig.Units = old;
            x = padL; y = padB; w = fp(3)-padL-padR; h = fp(4)-padT-padB;
            pos = [x/fp(3) y/fp(4) w/fp(3) h/fp(4)];
        end

        function pos = localSliderPos(self, fig, sliderW, padR, padT, padB)
            old = fig.Units; fig.Units = 'pixels';
            fp = fig.Position; fig.Units = old;
            x = fp(3)-padR-sliderW; y = padB;
            w = sliderW; h = fp(4)-padT-padB;
            pos = [x y w h];
        end

        function localOnHelpResize(self, padL,padR,padT,padB,sliderW,btn,contentH)
            if ~isvalid(self.helpFig), return; end
            % reposition axes & slider
            set(self.helpAx,'Position', self.localAxPos(self.helpFig, padL, padR+sliderW, padT, padB));
            set(self.helpSlider,'Position', self.localSliderPos(self.helpFig, sliderW, padR, padT, padB));
            % keep button bottom-right
            fp = getpixelposition(self.helpFig);
            set(btn,'Position',[fp(3)-90, 10, 80, 28]);
            % guard YLim within content
            cur = self.helpSlider.Value;
            maxScroll = max(0, contentH - 1);
            cur = min(max(cur,0), maxScroll);
            set(self.helpSlider,'Max',maxScroll,'Value',cur);
            set(self.helpAx,'YLim',[cur cur+1]);
        end

        function onHelpResize(self)
            if isempty(self.helpFig) || ~isvalid(self.helpFig), return; end
            h = self.helpFig;
            % Find the Close button (last created uicontrol)
            btns = findall(h, 'Style','pushbutton');
            if isempty(btns), return; end
            btn = btns(1);
            set(btn, 'Position', [h.Position(3)-90, 10, 80, 28]);
        end

        function placeHelpButton(self)
            if isempty(self.helpBtn) || ~ishandle(self.helpBtn), return; end
            oldUnits = get(self.f,'Units');
            set(self.f,'Units','pixels');
            fp = get(self.f,'Position'); % [x y w h]
            set(self.f,'Units',oldUnits);

            btnW = 26; btnH = 22;
            x = fp(1) - 2* btnW;
            y = fp(4) - 2* btnH;
            set(self.helpBtn, 'Units','pixels', 'Position', [x y btnW btnH]);
        end

        function call = generateVirtualBatCall(self, f_hi, f_lo, dur_s, fs, tailPct)
            % GENERATEVIRTUALBATCALL  Down-swept FM (f_hi -> f_lo) with windowing and padding.
            %
            % Inputs
            %   f_hi    : start frequency (Hz), typically higher (e.g., 90e3)
            %   f_lo    : end frequency (Hz), typically lower  (e.g., 25e3)
            %   dur_s   : duration in seconds (e.g., 0.005 for 5 ms)
            %   fs      : sampling rate (Hz). Pass [], 0, or omit to use self.fs
            %   tailPct : percent of dur_s to pad with zeros at head and tail (e.g., 100)
            %
            % Output
            %   call    : column vector, zero-padded at both ends

            if nargin < 5 || isempty(tailPct), tailPct = 0; end
            if nargin < 4 || isempty(fs) || fs <= 0, fs = self.fs; end

            % time vector (seconds)
            N  = max(1, round(dur_s * fs));
            t  = (0:N-1) / fs;

            % FM sweep: directly from f_hi to f_lo (quadratic to mimic bat calls)
            sig = chirp(t, f_lo, dur_s, f_hi, 'quadratic');
            sig = fliplr(sig);

            % Taper to reduce transients
            w = hann(N).';         % row window
            sig = sig .* w;

            % --- Optional emphasis filter (disabled by default) ---
            % fmax = mean([f_hi, f_lo]) - f_hi/3;            % your original heuristic
            % fb   = [0 2*[f_lo fmax f_hi] fs] / fs;         % normalized breakpoints
            % m    = [0 0 1 0 0];
            % [yb, ya] = yulewalk(4, fb, m);
            % sig = filtfilt(yb, ya, sig);

            % Normalize (optional; comment out if you prefer raw amplitude)
            % sig = sig ./ max(1e-12, max(abs(sig)));

            % Head/tail padding
            padN = round(N * tailPct / 100);
            call = [zeros(padN,1); sig(:); zeros(padN,1)];
        end
    end
    methods (Access = public)
        function exportAxestoPDF(self, ax, filename, varargin)
            % Export one axes exactly as shown (title, labels, ticks, legend/colorbar),
            % matching the GUI appearance. Defaults to PDF @ 300 dpi (image content).
            %
            % Usage:
            %   self.exportAxestoPDF(self.ax1, 'spectrogram.pdf');             % 300 dpi image-PDF
            %   self.exportAxestoPDF(self.ax2, 'diagram.pdf','Resolution',600) % 600 dpi
            %   self.exportAxestoPDF(self.ax1b,'wave.png');                    % PNG

            % --- on-screen size of the axes (pixels)
            axPix = getpixelposition(ax, true);  % [x y w h]
            w = max(1, round(axPix(3)));
            h = max(1, round(axPix(4)));

            % --- clean, invisible figure (same size), force OpenGL look
            tf = figure('Visible','off','Color','w','Units','pixels', ...
                'Position',[100 100 w h], 'MenuBar','none','ToolBar','none', ...
                'Renderer','opengl');

            % --- find linked legend/colorbar without creating them
            srcFig = ancestor(ax,'figure');
            lg = [];
            allLeg = findobj(srcFig,'Type','Legend');
            for k = 1:numel(allLeg)
                try
                    if isequal(allLeg(k).Axes, ax), lg = allLeg(k); break; end
                end
            end
            cb = findobj(srcFig,'Type','Colorbar','-and','Axes',ax);

            % --- copy axes (+legend together if it exists)
            if ~isempty(lg) && isvalid(lg)
                objs   = copyobj([ax lg], tf);        % keep legend tied to axes
                axCopy = objs(1);
            else
                axCopy = copyobj(ax, tf);
            end
            if ~isempty(cb) && isvalid(cb)
                copyobj(cb, tf);                       % copy colorbar if present
            end

            % --- preserve margins for title/labels/ticks using TightInset
            set(axCopy,'Units','normalized');
            ti  = get(axCopy,'TightInset');
            pos = [ti(1)+0.01, ti(2)+0.1, ...
                1- ti(3), ...
                1- ti(4)];
            set(axCopy,'Position', pos, 'ActivePositionProperty','position');
            % do NOT zero LooseInset — that would clip labels/titles

            % --- keep visual properties identical
            colormap(axCopy, colormap(ax));
            set(axCopy,'CLim',get(ax,'CLim'),'YDir',get(ax,'YDir'));
            drawnow;   % ensure fully rendered before export

            % --- choose export mode; default to image content @ 300 dpi for fidelity
            [~,~,ext] = fileparts(char(filename));
            hasRes = any(strcmpi(varargin(1:2:end),'Resolution'));
            hasCT  = any(strcmpi(varargin(1:2:end),'ContentType'));
            if strcmpi(ext,'.pdf')
                if ~hasCT,  varargin = [{'ContentType','image'} varargin]; end
                if ~hasRes, varargin = [{'Resolution',300} varargin];      end
            else
                if ~hasRes, varargin = [{'Resolution',300} varargin];      end
            end

            exportgraphics(tf, filename, varargin{:});
            close(tf);
        end

        function exportAxesExact(self, ax, filename, varargin)
            % Export one axes EXACTLY as seen on screen:
            % - preserves ticks, labels, titles, legends, colorbars
            % - preserves axes colormap and limits
            % - preserves figure background colour
            % - includes margins so labels/ticks are not cropped

            % --- Freeze the on-screen state
            srcFig = ancestor(ax,'figure');
            ax.Units = 'pixels';
            axOuter  = ax.OuterPosition;
            CMapAx   = colormap(ax);
            climAx   = get(ax,'CLim');
            alimAx   = get(ax,'ALim');
            bgColor  = get(srcFig,'Color');
            rendSrc  = get(srcFig,'Renderer');

            % --- New temp figure
            figTmp = figure('Visible','off','Units','pixels', ...
                'Position',[100 100 axOuter(3) axOuter(4)], ...
                'Color', bgColor,'MenuBar','none','ToolBar','none', ...
                'InvertHardcopy','off','Renderer',rendSrc);

            % --- Copy axes
            axCopy = copyobj(ax, figTmp);
            set(axCopy,'Units','normalized');

            % Reserve padding using TightInset
            ti  = get(axCopy,'TightInset');   % [L B R T]
            margin = 0.05; % 5% extra padding
            pos = [ti(1)+0.5*margin, ti(2)+1.6*margin, ...
                1 - ti(1) - ti(3) - 2*margin, ...
                1 - ti(2) - ti(4) - 2*margin];
            set(axCopy,'Position',pos);

            % --- Apply same colormap and limits
            colormap(axCopy,CMapAx);
            set(axCopy,'CLim',climAx);
            if ~isempty(alimAx), set(axCopy,'ALim',alimAx); end

            % --- Copy legend & colorbar
            lg = findobj(srcFig,'Type','Legend','-and','Axes',ax);
            if ~isempty(lg)&&isvalid(lg), copyobj(lg,figTmp); end
            cb = findobj(srcFig,'Type','Colorbar','-and','Axes',ax);
            if ~isempty(cb)&&isvalid(cb), copyobj(cb,figTmp); end

            drawnow;

            % --- Export
            [~,~,ext] = fileparts(char(filename));
            args = varargin;
            if strcmpi(ext,'.pdf')
                args = [{'ContentType','image','BackgroundColor','current'} args];
            else
                if ~any(strcmpi(args(1:2:end),'Resolution'))
                    args = [{'Resolution',300,'BackgroundColor','current'} args];
                else
                    args = [{'BackgroundColor','current'} args];
                end
            end

            exportgraphics(figTmp, filename, args{:});
            close(figTmp);
        end

        function setParams(self, params)
            % setParams - update GUI parameters programmatically
            %
            % Usage:
            %   h.setParams(struct('h',0.5,'theta_deg',30,'kr',6));

            % fetch all sliders (top->bottom order)
            vals = self.getSliders();
            sliders = findall(self.f, 'Style','slider');
            [~, idx] = sort(arrayfun(@(h) h.Position(2), sliders), 'descend');
            sliders = sliders(idx);

            % mapping of struct fields to slider order
            names = {'h','theta_deg','shore_dist','call_duration', ...
                'aperture','f_start_kHz','f_end_kHz','cutoff_db','kr'};

            for i = 1:numel(names)
                if isfield(params, names{i})
                    % convert units where needed
                    val = params.(names{i});
                    switch names{i}
                        case 'call_duration' % expects seconds
                            % leave as is
                        case {'f_start_kHz','f_end_kHz'} % expects kHz
                            % leave as is
                        otherwise
                            % all others are direct
                    end
                    set(sliders(i),'Value',val);
                end
            end

            % refresh once
            self.updatePlot();
        end
    end
end