%CONFIGURE_DYNAMIC_PLOTS Add time axes and desired-source motion to result figures.
%   figures: result handles; phiDS: [initial final] azimuths [rad].
%   SR: symbol rate [symbols/s]; N: symbol count.

function Configure_Dynamic_Plots(figures, phiDS, SR, N)
    duration = (N-1)/SR;
    for fig = figures(:).'
        if ~endsWith(fig.Name, ' | dynamic')
            fig.Name = [fig.Name ' | dynamic'];
        end
        resultLayout = findall(fig, 'Type', 'tiledlayout');
        resultLayout.Title.String = strrep(resultLayout.Title.String, ...
            'static scenario', 'dynamic scenario');
        if ~isappdata(fig, 'DynamicBaseSubtitle')
            setappdata(fig, 'DynamicBaseSubtitle', resultLayout.Subtitle.String);
        end
        baseSubtitle = getappdata(fig, 'DynamicBaseSubtitle');
        resultLayout.Subtitle.String = sprintf( ...
            '%s | DS: %.3g^{\\circ} to %.3g^{\\circ} | T = %.3g s', ...
            baseSubtitle, phiDS(1)*180/pi, phiDS(end)*180/pi, duration);

        for ax = findall(fig, 'Type', 'axes').'
            if any(strcmp(ax.XLabel.String, {'Iteration', 'Sample'}))
                sampleLimits = ax.XLim;
                for line = findall(ax, 'Type', 'line').'
                    line.XData = (line.XData-1)/SR;
                end
                for guide = findall(ax, 'Type', 'constantline').'
                    if strcmp(guide.InterceptAxis, 'x')
                        guide.Value = (guide.Value-1)/SR;
                    end
                end
                ax.XLim = (sampleLimits-1)/SR;
                ax.XTickMode = 'auto';
                ax.XLabel.String = 'Time [s]';
            end
            ax.Title.String = strrep(ax.Title.String, ...
                'dashed: Wiener', 'dashed: fixed Wiener');
            if startsWith(fig.Name, 'Wiener') && contains(ax.Title.String, 'constellation')
                ax.Title.String = '(c) Constellation | full run';
            end
        end

        patternAxes = findall(fig, 'Type', 'polaraxes');
        Add_Desired_Source_Motion(patternAxes, phiDS);
        if startsWith(fig.Name, 'Wiener')
            patternLabel = 'Pattern | fixed weights';
        else
            patternLabel = 'Pattern | final weights';
        end
        patternAxes.Title.String = strrep(patternAxes.Title.String, ...
            'Normalized radiation', patternLabel);
    end
end

function Add_Desired_Source_Motion(ax, phiDS)
    polarLegend = legend(ax);
    polarLegend.AutoUpdate = 'off';
    delete(findall(ax, '-regexp', 'Tag', 'DesiredDirection$|^DynamicDS'));
    patternLines = flipud(findall(ax, 'Type', 'line', '-regexp', 'Tag', 'Pattern$'));
    legendHandles = patternLines(:);
    legendLabels = arrayfun(@(line) line.DisplayName, patternLines(:), 'UniformOutput', false);
    sourceColor = [0.12 0.60 0.32];
    hold(ax, 'on');

    arcSize = max(2, min(2000, ceil(abs(diff(phiDS))*180/pi)+1));
    arcAngles = linspace(phiDS(1), phiDS(end), arcSize);
    pathRadius = ax.RLim(2) - 0.02*diff(ax.RLim);
    pathLine = polarplot(ax, arcAngles, pathRadius*ones(size(arcAngles)), ...
        'Color', sourceColor, 'LineWidth', 2, 'HandleVisibility', 'off', ...
        'Tag', 'DynamicDSPath');
    startMarker = polarplot(ax, phiDS(1), ax.RLim(2), 'o', ...
        'MarkerEdgeColor', sourceColor, 'MarkerFaceColor', 'w', ...
        'MarkerSize', 8, 'LineWidth', 1.6, 'LineStyle', 'none', 'Tag', 'DynamicDSStart');
    endMarker = polarplot(ax, phiDS(end), ax.RLim(2), '^', ...
        'MarkerEdgeColor', 'w', 'MarkerFaceColor', sourceColor, ...
        'MarkerSize', 9, 'LineStyle', 'none', 'Tag', 'DynamicDSEnd');
    legendHandles = [legendHandles; startMarker; endMarker];
    legendLabels{end+1, 1} = sprintf('DS start (%.3g^{\\circ})', phiDS(1)*180/pi);
    legendLabels{end+1, 1} = sprintf('DS end (%.3g^{\\circ})', phiDS(end)*180/pi);

    isMarkers = findall(ax, 'Type', 'line', '-regexp', 'Tag', 'InterferenceDirection$');
    if ~isempty(isMarkers)
        legendHandles = [legendHandles; isMarkers(:)];
        legendLabels = [legendLabels; repmat({'IS'}, numel(isMarkers), 1)];
        if isempty(findall(ax, 'Tag', 'InterferenceRay'))
            Plot_Interference_Rays(ax, isMarkers(1).ThetaData);
        end
    end
    legend(ax, legendHandles, legendLabels, 'Location', 'southoutside', ...
        'Orientation', 'horizontal', 'NumColumns', ceil(numel(legendHandles)/2), ...
        'FontSize', 9, 'Box', 'off', 'TextColor', [0.16 0.20 0.26]);
    limitListener = addlistener(ax, 'RLim', 'PostSet', ...
        @(~, ~) Update_Desired_Source_Motion(ax, pathLine, startMarker, endMarker));
    setappdata(ax, 'DynamicDSLimitListener', limitListener);
end

function Update_Desired_Source_Motion(ax, pathLine, startMarker, endMarker)
    if isgraphics(ax) && all(isgraphics([pathLine startMarker endMarker]))
        pathLine.RData(:) = ax.RLim(2) - 0.02*diff(ax.RLim);
        startMarker.RData = ax.RLim(2);
        endMarker.RData = ax.RLim(2);
    end
end
