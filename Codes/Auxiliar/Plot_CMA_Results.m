%PLOT_CMA_RESULTS Plot MSE, weights, constellation, and pattern.
%   y: output; eMSE: offline IQ error; W: K-by-N weight history.
%   constellation: reference symbols; phi: azimuth [rad].
%   patternCMA: normalized field; phiDS, phiIS: source azimuths [rad].
%   postStart: first constellation sample; mseWindow: averaging length.
%   algorithmName: figure label (default: CMA).
%   patternFloor_dB: radial display minimum (default: -35 dB).
%   SER: optional final symbol error rate after offline correction.
%   figCMA: figure handle.

function figCMA = Plot_CMA_Results(y, eMSE, W, constellation, phi, patternCMA, ...
    phiDS, phiIS, postStart, mseWindow, algorithmName, patternFloor_dB, SER)
    if nargin < 11
        algorithmName = 'CMA';
    end
    if nargin < 12
        patternFloor_dB = -35;
    end
    if nargin < 13
        SER = NaN;
    end
    validateattributes(patternFloor_dB, {'numeric'}, {'real', 'finite', 'scalar', '<', 0});
    plotTag = strrep(algorithmName, '-', '');

    N = numel(y);
    K = size(W, 1);
    plotIndex = unique(round(linspace(1, N, min(N, 10000))));
    constellationIndex = unique(round(linspace(postStart, N, ...
        min(8000, N-postStart+1))));
    assert(numel(eMSE) == N, 'MSE error must have one value per output sample.');
    mseEstimate = movmean(abs(eMSE(:)).^2, [mseWindow-1 0]);
    mse_dB = 10*log10(max(mseEstimate, realmin));

    gainColor = [0.08 0.38 0.68];
    referenceColor = [0.90 0.39 0.12];
    inkColor = [0.16 0.20 0.26];
    gridColor = [0.65 0.69 0.74];

    figCMA = figure('Name', [algorithmName ' Results'], 'NumberTitle', 'off', ...
        'Color', 'w', 'Position', [140 100 1250 900]);
    resultLayout = tiledlayout(figCMA, 2, 2, ...
        'TileSpacing', 'compact', 'Padding', 'loose');
    title(resultLayout, [algorithmName ' beamforming | static scenario'], ...
        'FontName', 'Arial', 'FontSize', 17, 'FontWeight', 'bold', ...
        'Color', inkColor);
    if ~isnan(SER)
        subtitle(resultLayout, sprintf('Final SER = %.3g', SER), ...
            'FontName', 'Arial', 'FontSize', 12, 'Color', inkColor);
    end

    errorAxes = nexttile(resultLayout, 1);
    plot(errorAxes, plotIndex, mse_dB(plotIndex), ...
        'Color', gainColor, 'LineWidth', 1.5, 'Tag', [plotTag 'MSE']);
    ylabel(errorAxes, 'MSE [dB]');
    title(errorAxes, sprintf('(a) MSE | %d-sample average', mseWindow), ...
        'FontSize', 12, 'Color', inkColor);
    xlim(errorAxes, [1 N]);
    xlabel(errorAxes, 'Iteration');

    weightAxes = nexttile(resultLayout, 2);
    weightLines = plot(weightAxes, plotIndex, abs(W(:,plotIndex)).', ...
        'LineWidth', 1.1);
    weightLabels = cell(K, 1);
    for k = 1:K
        weightLines(k).Tag = [plotTag 'Weight'];
        weightLabels{k} = sprintf('w_{%d}', k);
    end
    xlim(weightAxes, [1 N]);
    xlabel(weightAxes, 'Iteration');
    ylabel(weightAxes, '|w_k|');
    title(weightAxes, '(b) Coefficient magnitudes', ...
        'FontSize', 12, 'Color', inkColor);
    legend(weightAxes, weightLines, weightLabels, ...
        'Location', 'best', 'NumColumns', 2, 'Box', 'off', 'FontSize', 9, ...
        'TextColor', inkColor);

    constellationAxes = nexttile(resultLayout, 3);
    outputPoints = scatter(constellationAxes, real(y(constellationIndex)), ...
        imag(y(constellationIndex)), 12, gainColor, 'filled', ...
        'MarkerFaceAlpha', 0.22, 'MarkerEdgeAlpha', 0.22, ...
        'Tag', [plotTag 'Constellation']);
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
    title(constellationAxes, sprintf('(c) Constellation | last %d symbols', ...
        N-postStart+1), 'FontSize', 12, 'Color', inkColor);
    legend(constellationAxes, [outputPoints referencePoints], ...
        {[algorithmName ' output'], 'Reference'}, 'Location', 'southoutside', ...
        'Orientation', 'horizontal', 'Box', 'off', 'TextColor', inkColor);

    patternAxes = polaraxes(resultLayout);
    patternAxes.Layout.Tile = 4;
    pattern_dB = 20*log10(patternCMA);
    patternLine = polarplot(patternAxes, phi, pattern_dB, ...
        'Color', gainColor, 'LineWidth', 1.8, 'Tag', [plotTag 'Pattern']);
    Configure_Polar_dB_Axes(patternAxes, patternLine, [patternFloor_dB 0], 'dB');
    Plot_Interference_Rays(patternAxes, phiIS);
    hold(patternAxes, 'on');
    dsMarker = polarplot(patternAxes, phiDS, 0, 'o', ...
        'MarkerFaceColor', [0.12 0.60 0.32], 'MarkerEdgeColor', 'w', ...
        'MarkerSize', 8, 'LineStyle', 'none', 'Tag', [plotTag 'DesiredDirection']);
    legendHandles = [patternLine dsMarker];
    legendLabels = {'Normalized pattern', 'DS'};
    if ~isempty(phiIS)
        isMarkers = polarplot(patternAxes, phiIS, ...
            zeros(size(phiIS)), 'x', ...
            'Color', [0.80 0.16 0.20], 'LineWidth', 1.5, ...
            'MarkerSize', 8, 'LineStyle', 'none', 'Tag', [plotTag 'InterferenceDirection']);
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

    for ax = [errorAxes weightAxes constellationAxes]
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
