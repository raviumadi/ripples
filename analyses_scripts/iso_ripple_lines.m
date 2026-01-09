%% isoripple_contours_ripplestudio_geometry
% -------------------------------------------------------------------------
% PURPOSE
%   Compute and visualise a family of iso-ripple contours—curves of constant
%   spectral ripple spacing f_rpl—in the (beam angle, bat height) plane,
%   using the same geometric assumptions as the RippleStudio model.
%
%   Each contour represents all combinations of bat height (h) and beam
%   angle (θ) that yield the same two-path interference ripple frequency
%   at the receiver.
%
% -------------------------------------------------------------------------
% GEOMETRY (MATCHING RIPPLESTUDIO)
%   - Bat at height h above a flat reflecting surface.
%   - Emission axis pitched downward by θ relative to horizontal.
%   - Specular reflection occurs where the beam intersects the surface.
%
%   Horizontal run to reflection point:
%       X = h / tan(θ)
%
%   Reflected path:
%       Two identical legs of length
%           a = sqrt(h^2 + X^2)
%
%   Direct path:
%       Receiver lies at the same height as the bat; direct path is purely
%       horizontal:
%           L_dir = 2 X
%
%   Path-length difference:
%       ΔL = L_refl − L_dir = 2a − 2X
%
% -------------------------------------------------------------------------
% RIPPLE FREQUENCY
%   The spectral ripple spacing (null-to-null) produced by two-path
%   interference is:
%
%       f_rpl = c / ΔL
%
%   where c is the speed of sound. Small ΔL (shallow angles or low heights)
%   yield large f_rpl; large ΔL yield smaller ripple spacing.
%
%   A small epsilon is used in the denominator to guard against numerical
%   blow-ups as θ → 0.
%
% -------------------------------------------------------------------------
% PARAMETER SPACE
%   heights : 0.1 – 1.5 m (400 samples)
%   angles  : 5 – 60 deg (400 samples)
%
%   Meshgrids:
%     H    : bat height (m)
%     Adeg : beam angle (deg)
%
% -------------------------------------------------------------------------
% VISUALISATION
%   - Only contour lines are plotted (no background colour map).
%   - Each contour corresponds to a single ripple frequency f_rpl.
%   - Contours are labelled directly on the plot.
%
%   Iso-levels drawn (Hz):
%       1000, 2000, 3000, 4000, 5000
%
%   Axes:
%     x-axis : beam angle (degrees)
%     y-axis : bat height (metres)
%
%   Aspect ratio:
%     axis square (equal scaling of angle and height ranges)
%
% -------------------------------------------------------------------------
% OUTPUT
%   - Figure titled:
%       'Iso-$f_{ripple}$ Contours'
%   - Saved to:
%       ../manuscript/fig/isoripple_contour.(ext)
%     via saveFigure (helper function assumed on MATLAB path).
%
% -------------------------------------------------------------------------
% DEPENDENCIES
%   - formatLatex.m   (axes formatting helper)
%   - saveFigure.m   (figure export helper)
%
% -------------------------------------------------------------------------
% INTERPRETATION
%   - Lines indicate geometric degeneracy: many (h, θ) pairs produce the
%     same ripple spacing and are therefore acoustically indistinguishable
%     under a pure two-path interference model.
%   - Near-horizontal beams (small θ) cluster contours tightly, reflecting
%     the rapid increase of f_rpl as ΔL → 0.
%   - This plot provides an analytic map of the “ripple equivalence classes”
%     inherent to the RippleStudio geometry.
% -------------------------------------------------------------------------
clear; clc;

c = 343;                        % speed of sound (m/s)
heights = linspace(0.1, 1.5, 400);   % m
angles  = linspace(5, 60, 400);      % deg

[Adeg, H] = meshgrid(angles, heights);          % Adeg: deg, H: m
Th = deg2rad(Adeg);                              % radians
X  = H ./ tan(Th);                               % horizontal run to specular point
a  = hypot(H, X);                                % each leg of the reflected path
Ldir = 2 .* X;                                   % direct path (same y -> horizontal)
dL   = 2.*a - Ldir;                              % path difference
fRipple = c ./ max(dL, eps);                     % Hz (guard eps for tiny angles)

% Iso-levels to draw (Hz)
levels = 1000:1000:5000;

figure('Color','w', 'position', [100 100 600 600]); hold on; box on;
set(gcf,'Name','Iso-ripple contours (RippleStudio geometry)');

% Plot only lines (no background)
co = lines(numel(levels));
for i = 1:numel(levels)
    [C, hC] = contour(Adeg, H, fRipple, [levels(i) levels(i)], ...
        'LineWidth', 2, 'Color', co(i,:));
    % Optional: label each contour
    clabel(C, hC, 'Color', co(i,:), 'FontSize', 10, ...
        'LabelSpacing', 600, 'Interpreter', 'latex');
end

xlabel('Beam angle (deg)');
ylabel('Bat height (m)');
title('Iso-$f_{ripple}$ Contours');
legend(arrayfun(@(L) sprintf('$f_{rpl}$ %d Hz', L), levels, ...
       'UniformOutput', false), 'Interpreter', 'latex', 'Location','northeast');
formatLatex(gca);
axis square
xlim([min(angles) max(angles)]);
ylim([min(heights) max(heights)]);

%% Save Figure
fig_path = '../fig';
base_name = 'isoripple_contour';
saveFigure(gcf, fig_path, base_name);