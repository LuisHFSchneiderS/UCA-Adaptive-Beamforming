%PLOT_WIENER_RESULTS Plot Wiener error, weights, constellation, and pattern.
%   y, e: output/error; wo: fixed weights; Jmin: minimum training MSE.
%   constellation: reference symbols; phi: azimuth [rad].
%   patternWiener: normalized field; phiDS, phiIS: directions [rad].
%   mseWindow: averaging length; figWiener: figure handle.
%   patternFloor_dB: radial display minimum (default: -35 dB).
%   SER: optional final symbol error rate.

function figWiener = Plot_Wiener_Results(y, e, wo, Jmin, constellation, ...
    phi, patternWiener, phiDS, phiIS, mseWindow, patternFloor_dB, SER)
    if nargin < 11
        patternFloor_dB = -35;
    end
    if nargin < 12
        SER = NaN;
    end
    validateattributes(patternFloor_dB, {'numeric'}, {'real', 'finite', 'scalar', '<', 0});
    N = numel(y);
    K = numel(wo);
    plotIndex = unique(round(linspace(1, N, min(N, 10000))));
    constellationIndex = unique(round(linspace(1, N, min(N, 8000))));
    mseEstimate = movmean(abs(e).^2, [mseWindow-1 0]);
    mse_dB = 10*log10(max(mseEstimate, realmin));
    Jmin_dB = 10*log10(max(Jmin, realmin));

    gainColor = [0.08 0.38 0.68];
    referenceColor = [0.90 0.39 0.12];
    inkColor = [0.16 0.20 0.26];
    gridColor = [0.65 0.69 0.74];

    figWiener = figure('Name', 'Wiener Results', 'NumberTitle', 'off', ...
        'Color', 'w', 'Position', [110 80 1250 900]);
    resultLayout = tiledlayout(figWiener, 2, 2, ...
        'TileSpacing', 'compact', 'Padding', 'loose');
    title(resultLayout, 'Wiener beamforming | static scenario', ...
        'FontName', 'Arial', 'FontSize', 17, 'FontWeight', 'bold', ...
        'Color', inkColor);
    if ~isnan(SER)
        subtitle(resultLayout, sprintf('Final SER = %.3g', SER), ...
            'FontName', 'Arial', 'FontSize', 12, 'Color', inkColor);
    end

    mseAxes = nexttile(resultLayout, 1);
    mseLine = plot(mseAxes, plotIndex, mse_dB(plotIndex), ...
        'Color', gainColor, 'LineWidth', 1.5, 'Tag', 'WienerMSE');
    minimumLine = yline(mseAxes, Jmin_dB, '--', ...
        'Color', referenceColor, 'LineWidth', 1.3, 'Tag', 'WienerMinimumMSE');
    xlim(mseAxes, [1 N]);
    xlabel(mseAxes, 'Sample');
    ylabel(mseAxes, 'MSE [dB]');
    title(mseAxes, sprintf('(a) MSE estimate | %d-sample average', mseWindow), ...
        'FontSize', 12, 'Color', inkColor);
    legend(mseAxes, [mseLine minimumLine], {'Wiener output', 'J_{min}'}, ...
        'Location', 'best', 'Box', 'off', 'TextColor', inkColor);

    weightAxes = nexttile(resultLayout, 2);
    stem(weightAxes, 1:K, abs(wo), 'filled', 'Color', gainColor, ...
        'LineWidth', 1.5, 'Tag', 'WienerCoefficient');
    xlim(weightAxes, [0.5 K+0.5]);
    xticks(weightAxes, 1:K);
    xlabel(weightAxes, 'Coefficient index');
    ylabel(weightAxes, '|w_{o,k}|');
    title(weightAxes, '(b) Wiener coefficient magnitudes', ...
        'FontSize', 12, 'Color', inkColor);

    constellationAxes = nexttile(resultLayout, 3);
    outputPoints = scatter(constellationAxes, real(y(constellationIndex)), ...
        imag(y(constellationIndex)), 12, gainColor, 'filled', ...
        'MarkerFaceAlpha', 0.22, 'MarkerEdgeAlpha', 0.22, ...
        'Tag', 'WienerConstellation');
    hold(constellationAxes, 'on');
    referencePoints = scatter(constellationAxes, real(constellation), ...
        imag(constellation), 65, referenceColor, 'x', 'LineWidth', 1.6);
    constellationLimit = 1.2*max([abs(real(y(constellationIndex))); ...
        abs(imag(y(constellationIndex))); abs(real(constellation(:))); ...
        abs(imag(constellation(:))); 1]);
    axis(constellationAxes, 'equal');
    xlim(constellationAxes, [-constellationLimit constellationLimit]);
    ylim(constellationAxes, [-constellationLimit constellationLimit]);
    xlabel(constellationAxes, 'In-phase');
    ylabel(constellationAxes, 'Quadrature');
    title(constellationAxes, '(c) Wiener output constellation', ...
        'FontSize', 12, 'Color', inkColor);
    legend(constellationAxes, [outputPoints referencePoints], ...
        {'Wiener output', 'Reference'}, 'Location', 'southoutside', ...
        'Orientation', 'horizontal', 'Box', 'off', 'TextColor', inkColor);

    patternAxes = polaraxes(resultLayout);
    patternAxes.Layout.Tile = 4;
    pattern_dB = 20*log10(patternWiener);
    patternLine = polarplot(patternAxes, phi, pattern_dB, ...
        'Color', gainColor, 'LineWidth', 1.8, 'Tag', 'WienerPattern');
    Configure_Polar_dB_Axes(patternAxes, patternLine, [patternFloor_dB 0], 'dB');
    hold(patternAxes, 'on');
    dsMarker = polarplot(patternAxes, phiDS, 0, 'o', ...
        'MarkerFaceColor', [0.12 0.60 0.32], 'MarkerEdgeColor', 'w', ...
        'MarkerSize', 8, 'LineStyle', 'none', 'Tag', 'WienerDesiredDirection');
    legendHandles = [patternLine dsMarker];
    legendLabels = {'Normalized pattern', 'DS'};
    if ~isempty(phiIS)
        isMarkers = polarplot(patternAxes, phiIS, ...
            zeros(size(phiIS)), 'x', ...
            'Color', [0.80 0.16 0.20], 'LineWidth', 1.5, ...
            'MarkerSize', 8, 'LineStyle', 'none', ...
            'Tag', 'WienerInterferenceDirection');
        legendHandles(end+1) = isMarkers;
        legendLabels{end+1} = 'IS';
    end
    patternAxes.ThetaZeroLocation = 'right';
    patternAxes.ThetaDir = 'counterclockwise';
    patternAxes.ThetaTick = 0:30:330;
    patternAxes.RAxisLocation = 75;
    patternAxes.FontName = 'Arial';
    patternAxes.FontSize = 10;
    patternAxes.Color = 'w';
    patternAxes.RColor = inkColor;
    patternAxes.ThetaColor = inkColor;
    patternAxes.GridColor = gridColor;
    patternAxes.GridAlpha = 0.35;
    patternAxes.Clipping = 'on';
    grid(patternAxes, 'on');
    title(patternAxes, '(d) Normalized radiation | \theta = 90^{\circ}', ...
        'FontSize', 12, 'Color', inkColor);
    legend(patternAxes, legendHandles, legendLabels, ...
        'Location', 'southoutside', 'Orientation', 'horizontal', ...
        'Box', 'off', 'TextColor', inkColor);

    for ax = [mseAxes weightAxes constellationAxes]
        ax.FontName = 'Arial';
        ax.FontSize = 10.5;
        ax.Color = 'w';
        ax.XColor = inkColor;
        ax.YColor = inkColor;
        ax.GridColor = gridColor;
        ax.GridAlpha = 0.22;
        ax.Box = 'off';
        grid(ax, 'on');
    end
end
