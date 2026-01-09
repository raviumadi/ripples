%% interference_region_surface_plot.m
% -------------------------------------------------------------------------
% PURPOSE
%   Visualise the geometric “interference region” in (bat height, beam angle)
%   space for a two-path model, highlighting where the path-difference delay
%   is sufficiently small to produce closely spaced spectral ripples (or,
%   more generally, where a chosen delay criterion is met).
%
%   The script builds a 3D surface over:
%     - x-axis: bat height h (m)
%     - y-axis: beam angle θ (deg)
%     - z-axis: interference-point horizontal coordinate x_int (m)
%
%   A blue “base” surface shows the full geometry; a red overlay marks the
%   subset of parameter space satisfying the condition dt < 0.4 ms.
%
% -------------------------------------------------------------------------
% PARAMETERS
%   c          : speed of sound (m/s), default 343
%   shore_dist : microphone/shore horizontal distance (m), set to 10
%                (currently not used in the active mask; see below)
%
%   heights         : 0.1:0.01:1.5 m     (grid for h)
%   beam_angles_deg : 1:1:89 deg        (grid for θ)
%
%   Note:
%     θ is defined as a downward beam angle (relative to horizontal),
%     consistent with earlier RippleStudio-style geometry.
%
% -------------------------------------------------------------------------
% GEOMETRY (RIPPLESTUDIO-STYLE CONSTRUCTION)
%   Using a specular reflection point for a ray emitted at angle θ:
%
%     X_r   = h / tan(θ)          (horizontal distance from bat to reflection point)
%     x_int = 2 X_r               (horizontal coordinate of the “interference point”)
%
%   Reflected path:
%     a = h / sin(θ)              (bat → reflection point)
%     b = a                       (symmetric return to height h)
%     L_refl = a + b = 2a
%
%   Direct path (straight line to interference point):
%     L_dir = sqrt(x_int^2 + h^2)
%
%   Path difference and delay:
%     ΔL = L_refl − L_dir
%     dt = (ΔL / c) * 1000        (ms)
%
% -------------------------------------------------------------------------
% INTERFERENCE REGION MASK
%   valid_region is a Boolean mask defining the “ripple-capable” subset:
%
%     valid_region = (dt < 0.4)
%
%   Interpretation:
%     This threshold selects geometries where the two arrivals are separated
%     by less than 0.4 ms. Depending on your signal duration and bandwidth,
%     such small delays typically yield higher ripple spacing (in frequency)
%     and strong overlap in time.
%
%   Optional constraint (commented in the script):
%     valid_region = (dt < 0.4) & (x_int < shore_dist)
%   which additionally restricts the interference point to be “in front of”
%   the microphone/shore distance.
%
% -------------------------------------------------------------------------
% VISUALISATION
%   - Base surface (blue): x_int(h, θ) over the full grid.
%   - Overlay surface (red): same surface but only where valid_region is true
%     (values outside region set to NaN so they do not render).
%
%   Plot settings:
%     view(135, 30)          : angled 3D view
%     ylim([0 90])           : θ range
%     axis square            : square axes box
%     legend labels          : “No Ripples” (blue) vs “Ripples” (red overlay)
%
%  
%
% -------------------------------------------------------------------------
% OUTPUT
%   Saves figure using saveFigure helper:
%       saveFigure(gcf, fig_path, 'interference_regions')
%   to:
%       ../fig/interference_regions.(ext)
%
% -------------------------------------------------------------------------
% DEPENDENCIES
%   - saveFigure.m (custom helper)
%
% -------------------------------------------------------------------------
% NOTES / CAVEATS
%   - The “interference point” x_int here is a geometric construction tied
%     to the symmetric specular model; if you later define the receiver at a
%     fixed coordinate (e.g., mic at x = shore_dist, y = h_mic), recompute
%     L_dir and the feasibility region from explicit coordinates.
%   - The dt-threshold is an analysis choice. If you want the threshold to
%     correspond to a specific ripple spacing, you can convert:
%         f_rpl ≈ 1/dt   (Hz)  (for null-to-null spacing under a two-path model)
%     and threshold on f_rpl instead.
% -------------------------------------------------------------------------

% Parameters
fig_path = '../fig';
c = 343;                    % Speed of sound (m/s)
shore_dist = 10;             % Microphone horizontal distance
heights = 0.1:0.01:1.5;     % Bat heights (m)
beam_angles_deg = 1:1:89;   % Beam angles (deg)

[H, ThetaDeg] = meshgrid(heights, beam_angles_deg);
ThetaRad = deg2rad(ThetaDeg);

% Geometry
X_r = H ./ tan(ThetaRad);       % Horizontal distance to reflection
X_int = 2 * X_r;                % Total to interference point
a = H ./ sin(ThetaRad);         % First segment of reflection
b = a;                          % Symmetric reflection
total_path = a + b;
direct_path = sqrt(X_int.^2 + H.^2);
dd = total_path - direct_path;
dt = dd / c * 1000;             % Time delay in ms

% --- Condition mask
% valid_region = (dt < 0.4) & (X_int < shore_dist);
valid_region = (dt < 0.4);

% Create figure
figure('Position', [100 100 600 600]);
hold on;

% Base blue surface
base = surf(H, ThetaDeg, X_int, ...
    'EdgeColor', 'none', ...   % Gray grid lines
    'EdgeAlpha', 0.3, ...
    'FaceColor', [0 0 102]./255, ...      % BLUE
    'FaceAlpha', 1);

% Red overlay surface: set NaNs elsewhere
Z_red = X_int;
Z_red(~valid_region) = NaN;

overlay = surf(H, ThetaDeg, Z_red, ...
    'EdgeColor', 'none', ...   % Gray grid lines
    'EdgeAlpha', 0.3, ...
    'FaceColor', [153 0 0]./255, ...      % RED
    'FaceAlpha', 1);

% Plot formatting
xlabel('Bat Height h (m)', 'Fontsize', 14, 'Interpreter', 'latex');
ylabel('Beam Angle $\theta$', 'Fontsize', 14, 'Interpreter', 'latex');
zlabel('Interference Point $x_{int}$ (m)','Fontsize', 14, 'Interpreter', 'latex');
title('Interference Region', 'Fontsize', 18, 'Interpreter', 'latex');
view(135, 30);
ylim([0 90]);
grid on;
axis square
legend({'No Ripples', 'Ripples - Path Diff. $<$0.4ms'}, 'Location', 'east', 'Fontsize', 14, 'Interpreter', 'latex');
saveFigure(gcf, fig_path, 'interference_regions')
