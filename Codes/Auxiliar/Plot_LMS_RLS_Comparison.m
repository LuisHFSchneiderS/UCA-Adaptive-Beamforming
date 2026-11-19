%PLOT_LMS_RLS_COMPARISON Compare LMS and RLS beamforming results.
%   yLMS, yRLS: output IQ; eLMS, eRLS: reference errors.
%   constellation: reference symbols.
%   phi: azimuth [rad]; patternLMS, patternRLS: normalized fields.
%   phiDS, phiIS: source directions [rad]; postStart: first constellation sample.
%   mseWindow: averaging length; patternFloor_dB: radial minimum [dB].
%   SER: symbol error rates in [LMS, RLS] order.

function figComparison = Plot_LMS_RLS_Comparison(yLMS, yRLS, eLMS, eRLS, ...
    constellation, phi, patternLMS, patternRLS, ...
    phiDS, phiIS, postStart, mseWindow, patternFloor_dB, SER)
    if nargin < 13
        patternFloor_dB = -35;
    end
    if nargin < 14
        SER = [NaN NaN];
    end
    validateattributes(patternFloor_dB, {'numeric'}, {'real', 'finite', 'scalar', '<', 0});
    N = numel(yLMS);
    plotIndex = unique(round(linspace(1, N, min(N, 10000))));
    constellationIndex = unique(round(linspace(postStart, N, min(4000, N-postStart+1))));
    mse = movmean(abs([eLMS(:) eRLS(:)]).^2, [mseWindow-1 0], 1);
    mse_dB = 10*log10(max(mse, realmin));

    lmsColor = [0.08 0.38 0.68];
    rlsColor = [0.90 0.39 0.12];
    inkColor = [0.16 0.20 0.26];
    gridColor = [0.65 0.69 0.74];

    figComparison = figure('Name', 'LMS vs RLS', 'NumberTitle', 'off', ...
        'Color', 'w', 'Position', [80 60 1400 950]);
    resultLayout = tiledlayout(figComparison, 2, 2, ...
        'TileSpacing', 'compact', 'Padding', 'loose');
    title(resultLayout, 'LMS vs RLS | static scenario', ...
        'FontName', 'Arial', 'FontSize', 17, 'FontWeight', 'bold', 'Color', inkColor);
    if all(~isnan(SER))
        subtitle(resultLayout, sprintf('Final SER | LMS = %.3g | RLS = %.3g', SER), ...
            'FontName', 'Arial', 'FontSize', 12, 'Color', inkColor);
    end

    errorAxes = nexttile(resultLayout, 1, [1 2]);
    lmsMSE = plot(errorAxes, plotIndex, mse_dB(plotIndex, 1), ...
        'Color', lmsColor, 'LineWidth', 1.5, 'Tag', 'LMSComparisonMSE');
    hold(errorAxes, 'on');
    rlsMSE = plot(errorAxes, plotIndex, mse_dB(plotIndex, 2), '--', ...
        'Color', rlsColor, 'LineWidth', 1.5, 'Tag', 'RLSComparisonMSE');
    xlim(errorAxes, [1 N]);
    xlabel(errorAxes, 'Iteration'); ylabel(errorAxes, 'MSE [dB]');
    title(errorAxes, sprintf('(a) MSE | %d-sample average', mseWindow), ...
        'FontSize', 12, 'Color', inkColor);
    legend(errorAxes, [lmsMSE rlsMSE], {'LMS', 'RLS'}, ...
        'Location', 'best', 'Box', 'off', 'TextColor', inkColor);

    constellationAxes = nexttile(resultLayout, 3);
    lmsPoints = scatter(constellationAxes, real(yLMS(constellationIndex)), ...
        imag(yLMS(constellationIndex)), 12, lmsColor, 'filled', ...
        'MarkerFaceAlpha', 0.18, 'MarkerEdgeAlpha', 0.12, 'Tag', 'LMSComparisonConstellation');
    hold(constellationAxes, 'on');
    rlsPoints = scatter(constellationAxes, real(yRLS(constellationIndex)), ...
        imag(yRLS(constellationIndex)), 16, rlsColor, '+', ...
        'MarkerEdgeAlpha', 0.28, 'Tag', 'RLSComparisonConstellation');
    referencePoints = scatter(constellationAxes, real(constellation), imag(constellation), ...
        60, inkColor, 'x', 'LineWidth', 1.5);
    constellationLimit = 1.2*max([abs(real(yLMS(constellationIndex))); ...
        abs(imag(yLMS(constellationIndex))); abs(real(yRLS(constellationIndex))); ...
        abs(imag(yRLS(constellationIndex))); abs(real(constellation(:))); ...
        abs(imag(constellation(:))); 1]);
    axis(constellationAxes, 'equal');
    xlim(constellationAxes, [-constellationLimit constellationLimit]);
    ylim(constellationAxes, [-constellationLimit constellationLimit]);
    xlabel(constellationAxes, 'In-phase'); ylabel(constellationAxes, 'Quadrature');
    title(constellationAxes, sprintf('(b) Constellations | last %d symbols', N-postStart+1), ...
        'FontSize', 12, 'Color', inkColor);
    legend(constellationAxes, [lmsPoints rlsPoints referencePoints], ...
        {'LMS', 'RLS', 'Reference'}, 'Location', 'southoutside', ...
        'Orientation', 'horizontal', 'Box', 'off', 'TextColor', inkColor);

    patternAxes = polaraxes(resultLayout);
    patternAxes.Layout.Tile = 4;
    lmsPattern = polarplot(patternAxes, phi, 20*log10(patternLMS), ...
        'Color', lmsColor, 'LineWidth', 1.8, 'Tag', 'LMSComparisonPattern');
    hold(patternAxes, 'on');
    rlsPattern = polarplot(patternAxes, phi, 20*log10(patternRLS), '--', ...
        'Color', rlsColor, 'LineWidth', 1.8, 'Tag', 'RLSComparisonPattern');
    Configure_Polar_dB_Axes(patternAxes, lmsPattern, [patternFloor_dB 0], 'dB');
    Configure_Polar_dB_Axes(patternAxes, rlsPattern, [patternFloor_dB 0], 'dB');
    Plot_Interference_Rays(patternAxes, phiIS);
    dsMarker = polarplot(patternAxes, phiDS, 0, 'o', ...
        'MarkerFaceColor', [0.12 0.60 0.32], 'MarkerEdgeColor', 'w', ...
        'MarkerSize', 8, 'LineStyle', 'none', 'Tag', 'ComparisonDesiredDirection');
    legendHandles = [lmsPattern rlsPattern dsMarker];
    legendLabels = {'LMS', 'RLS', 'DS'};
    if ~isempty(phiIS)
        isMarkers = polarplot(patternAxes, phiIS, zeros(size(phiIS)), 'x', ...
            'Color', [0.80 0.16 0.20], 'LineWidth', 1.5, 'MarkerSize', 8, ...
            'LineStyle', 'none', 'Tag', 'ComparisonInterferenceDirection');
        legendHandles(end+1) = isMarkers;
        legendLabels{end+1} = 'IS';
    end
    patternAxes.ThetaZeroLocation = 'right';
    patternAxes.ThetaDir = 'counterclockwise';
    patternAxes.ThetaTick = 0:30:330;
    patternAxes.RAxisLocation = 75;
    patternAxes.FontName = 'Arial'; patternAxes.FontSize = 10;
    patternAxes.Color = 'w'; patternAxes.RColor = inkColor; patternAxes.ThetaColor = inkColor;
    patternAxes.GridColor = gridColor; patternAxes.GridAlpha = 0.35;
    grid(patternAxes, 'on');
    title(patternAxes, '(c) Normalized radiation | \theta = 90^{\circ}', ...
        'FontSize', 12, 'Color', inkColor);
    legend(patternAxes, legendHandles, legendLabels, 'Location', 'southoutside', ...
        'Orientation', 'horizontal', 'Box', 'off', 'TextColor', inkColor);

    for ax = [errorAxes constellationAxes]
        ax.FontName = 'Arial'; ax.FontSize = 10.5; ax.Color = 'w';
        ax.XColor = inkColor; ax.YColor = inkColor;
        ax.GridColor = gridColor; ax.GridAlpha = 0.22; ax.Box = 'off';
        grid(ax, 'on');
    end
end
