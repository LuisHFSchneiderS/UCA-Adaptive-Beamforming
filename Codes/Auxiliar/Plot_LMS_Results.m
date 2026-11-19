%PLOT_LMS_RESULTS Plot LMS error, weights, constellation, and radiation pattern.
%   y, e: N-by-1 output/error; W: K-by-N weight history.
%   constellation: reference symbols; phi: azimuth [rad].
%   patternLMS: normalized field; phiDS, phiIS: source azimuths [rad].
%   postStart: first constellation sample; mseWindow: averaging length.
%   wo: Wiener weights on the same input scale as W.
%   patternFloor_dB: radial display minimum (default: -35 dB).
%   SER: optional final symbol error rate.
%   algorithmName: plot label (default: 'LMS').
%   figLMS: figure handle.

function figLMS = Plot_LMS_Results(y, e, W, constellation, phi, patternLMS, ...
    phiDS, phiIS, postStart, mseWindow, wo, patternFloor_dB, SER, algorithmName)
    if nargin < 12
        patternFloor_dB = -35;
    end
    if nargin < 13
        SER = NaN;
    end
    if nargin < 14
        algorithmName = 'LMS';
    end
    plotTag = strrep(algorithmName, '-', '');
    validateattributes(patternFloor_dB, {'numeric'}, {'real', 'finite', 'scalar', '<', 0});
    N = numel(y);
    K = size(W, 1);
    plotIndex = unique(round(linspace(1, N, min(N, 10000))));
    constellationIndex = unique(round(linspace(postStart, N, ...
        min(8000, N-postStart+1))));
    mseEstimate = movmean(abs(e).^2, [mseWindow-1 0]);
    mse_dB = 10*log10(max(mseEstimate, realmin));

    gainColor = [0.08 0.38 0.68];
    referenceColor = [0.90 0.39 0.12];
    inkColor = [0.16 0.20 0.26];
    gridColor = [0.65 0.69 0.74];

    figLMS = figure('Name', [algorithmName ' Results'], 'NumberTitle', 'off', ...
        'Color', 'w', 'Position', [80 60 1250 900]);
    resultLayout = tiledlayout(figLMS, 2, 2, ...
        'TileSpacing', 'compact', 'Padding', 'loose');
    title(resultLayout, [algorithmName ' beamforming | static scenario'], ...
        'FontName', 'Arial', 'FontSize', 17, 'FontWeight', 'bold', ...
        'Color', inkColor);
    if ~isnan(SER)
        subtitle(resultLayout, sprintf('Final SER = %.3g', SER), ...
            'FontName', 'Arial', 'FontSize', 12, 'Color', inkColor);
    end

    mseAxes = nexttile(resultLayout, 1);
    plot(mseAxes, plotIndex, mse_dB(plotIndex), ...
        'Color', gainColor, 'LineWidth', 1.5, 'Tag', [plotTag 'MSE']);
    xlim(mseAxes, [1 N]);
    xlabel(mseAxes, 'Iteration');
    ylabel(mseAxes, 'MSE [dB]');
    title(mseAxes, sprintf('(a) MSE estimate | %d-sample average', mseWindow), ...
        'FontSize', 12, 'Color', inkColor);

    weightAxes = nexttile(resultLayout, 2);
    weightLines = plot(weightAxes, plotIndex, abs(W(:,plotIndex)).', ...
        'LineWidth', 1.1);
    wienerLines = gobjects(K, 1);
    weightLabels = cell(2*K, 1);
    for k = 1:K
        weightLines(k).Tag = [plotTag 'Weight'];
        wienerLines(k) = yline(weightAxes, abs(wo(k)), '--', ...
            'Color', weightLines(k).Color, 'LineWidth', 1.8, ...
            'Tag', 'WienerWeight');
        weightLabels{k} = sprintf('w_{%d}', k);
        weightLabels{K+k} = sprintf('w_{o,%d}', k);
    end
    xlim(weightAxes, [1 N]);
    xlabel(weightAxes, 'Iteration');
    ylabel(weightAxes, '|w_k|');
    title(weightAxes, '(b) Coefficients | dashed: Wiener', ...
        'FontSize', 12, 'Color', inkColor);
    legend(weightAxes, [weightLines; wienerLines], weightLabels, ...
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
    pattern_dB = 20*log10(patternLMS);
    patternLine = polarplot(patternAxes, phi, pattern_dB, ...
        'Color', gainColor, 'LineWidth', 1.8, 'Tag', [plotTag 'Pattern']);
    Configure_Polar_dB_Axes(patternAxes, patternLine, [patternFloor_dB 0], 'dB');
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
