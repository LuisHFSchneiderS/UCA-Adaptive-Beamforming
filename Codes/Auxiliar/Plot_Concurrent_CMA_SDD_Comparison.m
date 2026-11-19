%PLOT_CONCURRENT_CMA_SDD_COMPARISON Compare Concurrent and CMA-SDD results.
%   y: corrected IQ; e: offline IQ errors; constellation: reference symbols.
%   phi, phiDS, phiIS: azimuths [rad]; patterns: normalized fields.
%   postStart: first constellation sample; mseWindow: averaging length.
%   patternFloor_dB: radial minimum [dB]; SER: [concurrent CMA-DD, CMA-SDD] rates.

function figComparison = Plot_Concurrent_CMA_SDD_Comparison(yConcurrent, yCMASDD, ...
    eConcurrent, eCMASDD, constellation, phi, patternConcurrent, patternCMASDD, ...
    phiDS, phiIS, postStart, mseWindow, patternFloor_dB, SER)
    if nargin < 13
        patternFloor_dB = -35;
    end
    if nargin < 14
        SER = [NaN NaN];
    end
    N = numel(yConcurrent);
    plotIndex = unique(round(linspace(1, N, min(N, 10000))));
    constellationIndex = unique(round(linspace(postStart, N, min(4000, N-postStart+1))));
    mse = movmean(abs([eConcurrent(:) eCMASDD(:)]).^2, [mseWindow-1 0], 1);
    mse_dB = 10*log10(max(mse, realmin));

    concurrentColor = [0.08 0.38 0.68];
    cmasddColor = [0.90 0.39 0.12];
    inkColor = [0.16 0.20 0.26];
    gridColor = [0.65 0.69 0.74];

    figComparison = figure('Name', 'Concurrent vs CMA-SDD', 'NumberTitle', 'off', ...
        'Color', 'w', 'Position', [80 60 1400 950]);
    resultLayout = tiledlayout(figComparison, 2, 2, ...
        'TileSpacing', 'compact', 'Padding', 'loose');
    title(resultLayout, 'Concurrent vs CMA-SDD | static scenario', ...
        'FontName', 'Arial', 'FontSize', 17, 'FontWeight', 'bold', 'Color', inkColor);
    if all(~isnan(SER))
        subtitle(resultLayout, sprintf('Final SER | Concurrent = %.3g | CMA-SDD = %.3g', SER), ...
            'FontName', 'Arial', 'FontSize', 12, 'Color', inkColor);
    end

    errorAxes = nexttile(resultLayout, 1, [1 2]);
    concurrentMSE = plot(errorAxes, plotIndex, mse_dB(plotIndex, 1), ...
        'Color', concurrentColor, 'LineWidth', 1.5, 'Tag', 'ConcurrentComparisonMSE');
    hold(errorAxes, 'on');
    cmasddMSE = plot(errorAxes, plotIndex, mse_dB(plotIndex, 2), '--', ...
        'Color', cmasddColor, 'LineWidth', 1.5, 'Tag', 'CMASDDComparisonMSE');
    xlim(errorAxes, [1 N]);
    xlabel(errorAxes, 'Iteration'); ylabel(errorAxes, 'MSE [dB]');
    title(errorAxes, sprintf('(a) MSE | %d-sample average', mseWindow), ...
        'FontSize', 12, 'Color', inkColor);
    legend(errorAxes, [concurrentMSE cmasddMSE], {'Concurrent', 'CMA-SDD'}, ...
        'Location', 'best', 'Box', 'off', 'TextColor', inkColor);

    constellationAxes = nexttile(resultLayout, 3);
    concurrentPoints = scatter(constellationAxes, real(yConcurrent(constellationIndex)), ...
        imag(yConcurrent(constellationIndex)), 12, concurrentColor, 'filled', ...
        'MarkerFaceAlpha', 0.18, 'MarkerEdgeAlpha', 0.12, 'Tag', 'ConcurrentComparisonConstellation');
    hold(constellationAxes, 'on');
    cmasddPoints = scatter(constellationAxes, real(yCMASDD(constellationIndex)), ...
        imag(yCMASDD(constellationIndex)), 16, cmasddColor, '+', ...
        'MarkerEdgeAlpha', 0.28, 'Tag', 'CMASDDComparisonConstellation');
    referencePoints = scatter(constellationAxes, real(constellation), imag(constellation), ...
        60, inkColor, 'x', 'LineWidth', 1.5);
    constellationLimit = 1.2*max([abs(real(yConcurrent(constellationIndex))); ...
        abs(imag(yConcurrent(constellationIndex))); abs(real(yCMASDD(constellationIndex))); ...
        abs(imag(yCMASDD(constellationIndex))); abs(real(constellation(:))); ...
        abs(imag(constellation(:))); 1]);
    axis(constellationAxes, 'equal');
    xlim(constellationAxes, [-constellationLimit constellationLimit]);
    ylim(constellationAxes, [-constellationLimit constellationLimit]);
    xlabel(constellationAxes, 'In-phase'); ylabel(constellationAxes, 'Quadrature');
    title(constellationAxes, sprintf('(b) Constellations | last %d symbols', N-postStart+1), ...
        'FontSize', 12, 'Color', inkColor);
    legend(constellationAxes, [concurrentPoints cmasddPoints referencePoints], ...
        {'Concurrent', 'CMA-SDD', 'Reference'}, 'Location', 'southoutside', ...
        'Orientation', 'horizontal', 'Box', 'off', 'TextColor', inkColor);

    patternAxes = polaraxes(resultLayout);
    patternAxes.Layout.Tile = 4;
    concurrentPattern = polarplot(patternAxes, phi, 20*log10(patternConcurrent), ...
        'Color', concurrentColor, 'LineWidth', 1.8, 'Tag', 'ConcurrentComparisonPattern');
    hold(patternAxes, 'on');
    cmasddPattern = polarplot(patternAxes, phi, 20*log10(patternCMASDD), '--', ...
        'Color', cmasddColor, 'LineWidth', 1.8, 'Tag', 'CMASDDComparisonPattern');
    Configure_Polar_dB_Axes(patternAxes, concurrentPattern, [patternFloor_dB 0], 'dB');
    Configure_Polar_dB_Axes(patternAxes, cmasddPattern, [patternFloor_dB 0], 'dB');
    Plot_Interference_Rays(patternAxes, phiIS);
    dsMarker = polarplot(patternAxes, phiDS, 0, 'o', ...
        'MarkerFaceColor', [0.12 0.60 0.32], 'MarkerEdgeColor', 'w', ...
        'MarkerSize', 8, 'LineStyle', 'none', 'Tag', 'ComparisonDesiredDirection');
    legendHandles = [concurrentPattern cmasddPattern dsMarker];
    legendLabels = {'Concurrent', 'CMA-SDD', 'DS'};
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
