%ANIMATE_BEAMFORMING_RESULTS Visual playback of existing base-workspace results.
%   method: Wiener, LMS, RLS, CMA, or Concurrent.
%   Duration: replay length [s]; FPS: refresh rate; WindowSize: SER/IQ window.
%   AutoPlay: start immediately. No filter or workspace data is changed.
%   GIFFile: optional GIF output path; Duration and FPS also control export.

function fig = Animate_Beamforming_Results(method, varargin)
    method = validatestring(method, {'Wiener', 'LMS', 'RLS', 'CMA', 'Concurrent'});
    parser = inputParser;
    addParameter(parser, 'Duration', 20, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
    addParameter(parser, 'FPS', 20, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && x <= 60);
    addParameter(parser, 'WindowSize', 2000, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && x == fix(x));
    addParameter(parser, 'AutoPlay', true, @(x) islogical(x) && isscalar(x));
    addParameter(parser, 'GIFFile', '', @(x) (ischar(x) && (isrow(x) || isempty(x))) || (isstring(x) && isscalar(x)));
    parse(parser, varargin{:});
    options = parser.Results;

    if any(strcmp(method, {'CMA', 'Concurrent'}))
        y = Read_Base(['y' method '_corrected']);
        e = Read_Base(['eMSE_' method]);
        decisions = Read_Base(['Dstream' method '_corrected']);
    else
        y = Read_Base(['y' method]);
        e = Read_Base(['e' method]);
        decisions = Read_Base(['Dstream' method]);
    end
    y = y(:); e = e(:); decisions = decisions(:);
    reference = Read_Base('Dstream');
    reference = reference(:);
    N = numel(y);
    assert(N >= 2 && numel(e) == N && numel(decisions) == N && numel(reference) == N, ...
        'Workspace outputs and symbol sequences must have matching lengths.');
    fixedWeights = strcmp(method, 'Wiener');
    if fixedWeights
        weights = Read_Base('woWiener');
        weights = weights(:);
    else
        weights = Read_Base(['W' method]);
        assert(size(weights, 2) == N, 'Weight history must have one column per sample.');
    end
    K = size(weights, 1);
    SR = Read_Base('SR');
    phi0 = Read_Base('phiDS_0');
    omega = Read_Base('ang_vel');
    constellation = Read_Base('constellation');
    validateattributes(SR, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
    validateattributes(omega, {'numeric'}, {'scalar', 'real', 'finite'});
    inputScale = Read_Base('inputScale');
    Z = Read_Base('Z'); ZL = Read_Base('ZL'); C = Read_Base('C');
    L = Read_Base('L'); wavelength = Read_Base('lambda'); r = Read_Base('r');
    phiCount = Read_Base('phi_numel');
    assert(isequal(size(Z), [K K]) && isequal(size(ZL), [K K]) && size(C,1) == K, ...
        'Array matrices and saved weights have incompatible dimensions.');
    mseWindow = min(N, Read_Base('mseWindow', 1000));
    radialFloor = Read_Base('patternFloor_dB', -30);
    phiIS = Read_Base('phiIS', []);
    phiIS = phiIS(1:Read_Base('Num_IS', numel(phiIS)));
    wienerWeights = Read_Base('woWiener', []);
    windowSize = min(N, options.WindowSize);
    phi = linspace(0, 2*pi, phiCount).';
    endPhi = phi0 + omega*(N-1)/SR;
    currentMap = ((Z + ZL)'\ZL')/(inputScale*sin(pi*L/wavelength));
    fieldMap = complex(zeros(phiCount, K));
    for k = 1:K
        fieldMap(:, k) = Calculate_Radiation_Pattern(currentMap(:, k), r, phiCount, L, C, wavelength);
    end
    finalField = fieldMap*weights(:, end);
    finalPeak = max(abs(finalField));
    if isfinite(finalPeak) && finalPeak > 0
        finalPattern = abs(finalField)/finalPeak;
    else
        finalPattern = ones(phiCount, 1);
    end
    tailStart = floor(0.8*N) + 1;
    validationMask = Read_Base('validationMask', ...
        Read_Base(['validationMask' method], true(N,1)));
    validationMask = validationMask(:);
    symbolErrors = decisions ~= reference;
    finalSER = mean(symbolErrors(validationMask));
    templateY = y;
    if any(~isfinite(templateY)), templateY(~isfinite(templateY)) = 0; end
    if fixedWeights
        fig = Plot_Wiener_Results(templateY, e, weights, Read_Base('JminWiener', mean(abs(e).^2)), ...
            constellation, phi, finalPattern, phi0, phiIS, mseWindow, radialFloor, finalSER);
    elseif any(strcmp(method, {'LMS', 'RLS'})) && numel(wienerWeights) == K
        fig = Plot_LMS_Results(templateY, e, weights, constellation, phi, finalPattern, ...
            phi0, phiIS, tailStart, mseWindow, wienerWeights, radialFloor, finalSER, method);
    else
        fig = Plot_CMA_Results(templateY, e, weights, constellation, phi, finalPattern, ...
            phi0, phiIS, tailStart, mseWindow, method, radialFloor, finalSER);
    end
    Configure_Dynamic_Plots(fig, [phi0 endPhi], SR, N);
    set(findall(fig, '-regexp', 'Tag', '^DynamicDS(Path|Start|End)$'), 'Visible', 'off');
    fig.Name = [method ' | dynamic replay'];
    layout = findall(fig, 'Type', 'tiledlayout');
    layout.Title.String = [method ' beamforming | dynamic replay'];
    layout.Units = 'normalized';
    layout.OuterPosition = [0 0.075 1 0.925];
    mseLine = findall(fig, 'Tag', [method 'MSE']);
    points = findall(fig, 'Tag', [method 'Constellation']);
    patternLine = findall(fig, 'Tag', [method 'Pattern']);
    polarAx = patternLine.Parent;
    mseAx = mseLine.Parent;
    constellationAx = points.Parent;
    plotIndex = unique(round(linspace(1, N, min(N, 10000)))).';
    mse = movmean(abs(e).^2, [mseWindow-1 0]);
    mseDB = 10*log10(max(mse, realmin));
    validationCount = movsum(double(validationMask), [windowSize-1 0]);
    serEstimate = movsum(double(symbolErrors & validationMask), [windowSize-1 0])./validationCount;
    serAx = nexttile(layout, 2);
    legend(serAx, 'off');
    cla(serAx);
    serLine = plot(serAx, NaN, NaN, 'Color', [0.08 0.38 0.68], ...
        'LineWidth', 1.5, 'Tag', [method 'SER']);
    hold(serAx, 'on');
    serZeroLine = plot(serAx, NaN, NaN, '.', 'Color', [0.08 0.38 0.68], ...
        'MarkerSize', 6, 'LineStyle', 'none', 'Tag', [method 'SERZero']);
    serResolution = 1/windowSize;
    zeroSERLevel = serResolution/10;
    serUpper = min(1, 10^ceil(log10(max([serResolution; serEstimate(isfinite(serEstimate))]))));
    serTicks = unique([zeroSERLevel serResolution ...
        10.^(ceil(log10(serResolution)):floor(log10(serUpper)))]);
    serTickLabels = arrayfun(@(value) sprintf('%.3g', value), serTicks, 'UniformOutput', false);
    serTickLabels{1} = '0 errors';
    xlim(serAx, [0 (N-1)/SR]);
    serAx.YScale = 'log';
    ylim(serAx, [0.8*zeroSERLevel 1.15*serUpper]);
    serAx.YTick = serTicks;
    serAx.YTickLabel = serTickLabels;
    serAx.XTickMode = 'auto';
    serAx.YMinorGrid = 'on';
    serAx.MinorGridAlpha = 0.08;
    serAx.FontName = 'Arial'; serAx.FontSize = 10.5;
    serAx.Color = 'w'; serAx.Box = 'off';
    serAx.XColor = [0.16 0.20 0.26]; serAx.YColor = serAx.XColor;
    serAx.GridColor = [0.65 0.69 0.74]; serAx.GridAlpha = 0.22;
    grid(serAx, 'on');
    xlabel(serAx, 'Time [s]'); ylabel(serAx, 'SER (log scale)');
    title(serAx, sprintf('(b) Validation SER | %d-symbol window', windowSize), ...
        'FontSize', 12, 'Color', [0.16 0.20 0.26]);
    animatedAxes = [mseAx serAx];
    cursors = gobjects(numel(animatedAxes), 1);
    for k = 1:numel(animatedAxes)
        guides = findall(animatedAxes(k), 'Type', 'constantline');
        for guide = guides(:).'
            if strcmp(guide.InterceptAxis, 'x'), delete(guide); end
        end
        cursors(k) = xline(animatedAxes(k), 0, ':', 'Color', [0.16 0.20 0.26], ...
            'LineWidth', 1.3, 'HandleVisibility', 'off', 'Tag', 'ReplayTimeCursor');
    end
    delete(getappdata(patternLine, 'RadialLimitListener'));
    patternListener = addlistener(polarAx, 'RLim', 'PostSet', @(~, ~) Refresh_Pattern());
    setappdata(patternLine, 'RadialLimitListener', patternListener);
    rawPattern = nan(phiCount, 1);
    sourceColor = [0.12 0.60 0.32];
    hold(polarAx, 'on');
    currentRay = polarplot(polarAx, [phi0 phi0], polarAx.RLim, ...
        'Color', sourceColor, 'LineWidth', 1.2, 'HandleVisibility', 'off', 'Tag', 'ReplayDSRay');
    currentRay.Color = [sourceColor 0.35];
    currentDS = polarplot(polarAx, phi0, polarAx.RLim(2), 'o', ...
        'MarkerFaceColor', sourceColor, 'MarkerEdgeColor', 'w', 'MarkerSize', 9, ...
        'LineStyle', 'none', 'Tag', 'ReplayDSCurrent');
    legendHandles = [patternLine; currentDS];
    legendLabels = {'Normalized pattern'; 'DS now'};
    isMarkers = findall(polarAx, 'Type', 'line', '-regexp', 'Tag', 'InterferenceDirection$');
    if ~isempty(isMarkers)
        legendHandles = [legendHandles; isMarkers(:)];
        legendLabels = [legendLabels; repmat({'IS'}, numel(isMarkers), 1)];
    end
    legend(polarAx, legendHandles, legendLabels, 'Location', 'southoutside', ...
        'Orientation', 'horizontal', 'NumColumns', 3, 'FontSize', 9, 'Box', 'off');

    playButton = uicontrol(fig, 'Style', 'togglebutton', 'Units', 'normalized', ...
        'Position', [0.02 0.017 0.07 0.035], 'String', 'Play', ...
        'Callback', @Toggle_Play, 'Tag', 'ReplayPlay');
    restartButton = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
        'Position', [0.10 0.017 0.075 0.035], 'String', 'Restart', ...
        'Callback', @(~, ~) Seek_Sample(1), 'Tag', 'ReplayRestart');
    saveButton = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
        'Position', [0.185 0.017 0.09 0.035], 'String', 'Save GIF', ...
        'Callback', @Choose_GIF, 'Tag', 'ReplaySaveGIF');
    slider = uicontrol(fig, 'Style', 'slider', 'Units', 'normalized', ...
        'Position', [0.30 0.022 0.37 0.025], 'Min', 0, 'Max', 1, 'Value', 0, ...
        'Callback', @(control, ~) Seek_Sample(1 + round(control.Value*(N-1))), 'Tag', 'ReplaySlider');
    speedControl = uicontrol(fig, 'Style', 'popupmenu', 'Units', 'normalized', ...
        'Position', [0.69 0.017 0.08 0.035], 'String', {'0.5x', '1x', '2x', '4x'}, ...
        'Value', 2, 'Callback', @Change_Speed, 'Tag', 'ReplaySpeed');
    counter = uicontrol(fig, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.79 0.017 0.19 0.035], 'BackgroundColor', 'w', ...
        'FontName', 'Arial', 'FontSize', 10, 'Tag', 'ReplayCounter');
    sample = 1;
    anchorFraction = 0;
    playbackClock = tic;
    isClosing = false;
    isExporting = false;
    cancelExport = false;
    playbackTimer = timer('ExecutionMode', 'fixedSpacing', 'Period', round(1000/options.FPS)/1000, ...
        'BusyMode', 'drop', 'Name', [method ' visual replay'], ...
        'TimerFcn', @Advance_Frame, 'ErrorFcn', @Playback_Error);
    setappdata(fig, 'ReplayTimer', playbackTimer);
    setappdata(fig, 'ReplaySeek', @Seek_Sample);
    setappdata(fig, 'ReplayMethod', method);
    setappdata(fig, 'ReplaySaveGIF', @Save_GIF);
    fig.CloseRequestFcn = @Close_Replay;
    Render_Frame();
    if ~isempty(char(options.GIFFile))
        Save_GIF(char(options.GIFFile));
    end
    if options.AutoPlay && isgraphics(fig)
        playButton.Value = 1;
        Toggle_Play();
    end

    function Render_Frame()
        if isClosing || ~isgraphics(fig), return; end
        t = (sample-1)/SR;
        angle = phi0 + omega*t;
        selected = plotIndex(plotIndex <= sample);
        if selected(end) ~= sample, selected(end+1) = sample; end
        set(mseLine, 'XData', (selected-1)/SR, 'YData', mseDB(selected));
        displayedSER = serEstimate(selected);
        zeroSER = displayedSER == 0;
        positiveSER = displayedSER;
        positiveSER(zeroSER) = NaN;
        set(serLine, 'XData', (selected-1)/SR, 'YData', positiveSER);
        setappdata(serLine, 'RawSER', displayedSER);
        serLine.DataTipTemplate.DataTipRows = [ ...
            dataTipTextRow('Time [s]', (selected-1)/SR); ...
            dataTipTextRow('SER', displayedSER)];
        zeroTimes = (selected(zeroSER)-1)/SR;
        set(serZeroLine, 'XData', zeroTimes, 'YData', zeroSERLevel*ones(size(zeroTimes)));
        serZeroLine.DataTipTemplate.DataTipRows = [ ...
            dataTipTextRow('Time [s]', zeroTimes); ...
            dataTipTextRow('SER', zeros(size(zeroTimes)))];
        for cursor = cursors(:).', cursor.Value = t; end
        recent = max(1, sample-windowSize+1):sample;
        shown = recent(unique(round(linspace(1, numel(recent), min(1500, numel(recent))))));
        displayedY = y(shown);
        displayedY = displayedY(isfinite(displayedY));
        set(points, 'XData', real(displayedY), 'YData', imag(displayedY));
        axisLimit = 1.2*max([abs(real(displayedY(:))); abs(imag(displayedY(:))); ...
            abs(real(constellation(:))); abs(imag(constellation(:))); 1]);
        xlim(constellationAx, [-axisLimit axisLimit]);
        ylim(constellationAx, [-axisLimit axisLimit]);
        constellationAx.Title.String = sprintf('(c) Constellation | last %d symbols', numel(recent));
        if fixedWeights, currentWeights = weights; else, currentWeights = weights(:, sample); end
        field = fieldMap*currentWeights;
        peak = max(abs(field));
        if all(isfinite(field)) && peak > 0
            rawPattern = 20*log10(max(abs(field)/peak, realmin));
            patternLine.Visible = 'on';
            if fixedWeights, status = 'fixed weights'; else, status = 'current weights'; end
        else
            rawPattern(:) = NaN;
            patternLine.Visible = 'off';
            if all(isfinite(field)), status = 'zero weights'; else, status = 'nonfinite weights'; end
        end
        polarAx.Title.String = ['(d) Pattern | ' status ' | \theta = 90^{\circ}'];
        Refresh_Pattern();
        currentDS.ThetaData = angle;
        currentDS.RData = polarAx.RLim(2);
        currentRay.ThetaData = [angle angle];
        currentRay.RData = polarAx.RLim;
        windowSER = serEstimate(sample);
        serLabel = sprintf('%.3g', windowSER);
        if isnan(windowSER), serLabel = 'n/a'; end
        windowMSE = mean(abs(e(recent)).^2);
        layout.Subtitle.String = sprintf( ...
            't = %.4g s | DS = %.3g^{\\circ} | Window SER = %s | MSE = %.3g dB', ...
            t, angle*180/pi, serLabel, 10*log10(max(windowMSE, realmin)));
        slider.Value = (sample-1)/(N-1);
        counter.String = sprintf('%d / %d', sample, N);
        setappdata(fig, 'ReplaySample', sample);
        setappdata(fig, 'ReplayTime', t);
        setappdata(fig, 'ReplayWindowSER', windowSER);
        drawnow limitrate;
    end

    function Refresh_Pattern()
        if isClosing || ~isgraphics(patternLine), return; end
        setappdata(patternLine, 'RawLevel_dB', rawPattern);
        patternLine.RData = min(max(rawPattern, polarAx.RLim(1)), polarAx.RLim(2));
        patternLine.DataTipTemplate.DataTipRows = [ ...
            dataTipTextRow('Azimuth [deg]', phi*180/pi); ...
            dataTipTextRow('Level [dB]', rawPattern)];
        tickStep = 5*max(1, ceil(diff(polarAx.RLim)/50));
        polarAx.RTick = unique([polarAx.RLim(1), ...
            ceil(polarAx.RLim(1)/tickStep)*tickStep:tickStep:polarAx.RLim(2), polarAx.RLim(2)]);
        currentDS.RData = polarAx.RLim(2);
        currentRay.RData = polarAx.RLim;
    end

    function Seek_Sample(index)
        sample = max(1, min(N, round(index)));
        anchorFraction = (sample-1)/(N-1);
        playbackClock = tic;
        Render_Frame();
    end

    function Toggle_Play(varargin)
        if playButton.Value
            if sample == N, Seek_Sample(1); end
            anchorFraction = (sample-1)/(N-1);
            playbackClock = tic;
            playButton.String = 'Pause';
            start(playbackTimer);
        else
            stop(playbackTimer);
            playButton.String = 'Play';
        end
    end

    function Change_Speed(varargin)
        anchorFraction = (sample-1)/(N-1);
        playbackClock = tic;
    end

    function Advance_Frame(varargin)
        if isClosing || ~isgraphics(fig), return; end
        speeds = [0.5 1 2 4];
        fraction = min(1, anchorFraction + toc(playbackClock)*speeds(speedControl.Value)/options.Duration);
        sample = 1 + round(fraction*(N-1));
        Render_Frame();
        if fraction >= 1 && ~isClosing && isgraphics(fig)
            stop(playbackTimer);
            playButton.Value = 0;
            playButton.String = 'Play';
        end
    end

    function Playback_Error(varargin)
        if ~isClosing && isgraphics(playButton)
            playButton.Value = 0;
            playButton.String = 'Play';
        end
    end

    function Choose_GIF(varargin)
        if isExporting
            cancelExport = true;
            return;
        end
        codeFolder = fileparts(fileparts(mfilename('fullpath')));
        defaultFolder = fullfile(fileparts(codeFolder), 'Figuras');
        if ~isfolder(defaultFolder), defaultFolder = pwd; end
        [file, folder] = uiputfile({'*.gif', 'Animated GIF (*.gif)'}, ...
            'Save animation as GIF', fullfile(defaultFolder, [method '_Dynamic.gif']));
        if isequal(file, 0), return; end
        try
            Save_GIF(fullfile(folder, file));
        catch exception
            if isgraphics(fig), errordlg(exception.message, 'GIF export'); end
        end
    end

    function Save_GIF(file)
        if isExporting || isClosing, return; end
        [folder, name, extension] = fileparts(char(file));
        assert(~isempty(name), 'Provide a GIF filename.');
        if isempty(extension), extension = '.gif'; end
        assert(strcmpi(extension, '.gif'), 'The output extension must be .gif.');
        if isempty(folder), folder = pwd; end
        if ~isfolder(folder)
            [success, message] = mkdir(folder);
            assert(success, '%s', message);
        end
        [~, attributes] = fileattrib(folder);
        file = fullfile(attributes.Name, [name extension]);
        frameFile = [tempname '.png'];
        temporaryGIF = [tempname(attributes.Name) '.gif'];
        savedSample = sample;
        wasPlaying = strcmp(playbackTimer.Running, 'on');
        savedResize = fig.Resize;
        stop(playbackTimer);
        playButton.Value = 0;
        playButton.String = 'Play';
        isExporting = true;
        cancelExport = false;
        controls = [playButton restartButton slider speedControl];
        set(controls, 'Enable', 'off');
        saveButton.String = 'Cancel GIF';
        fig.Resize = 'off';
        cleanup = onCleanup(@() Finish_Export(frameFile, temporaryGIF, ...
            savedSample, wasPlaying, savedResize, controls));
        frameCount = max(2, round(options.Duration*options.FPS));
        frameSamples = round(linspace(1, N, frameCount));
        centiseconds = max(frameCount, round(100*options.Duration));
        delays = diff(round(linspace(0, centiseconds, frameCount+1)))/100;
        figurePixels = getpixelposition(fig);
        imageHeight = round(1000*figurePixels(4)*layout.OuterPosition(4)/figurePixels(3));
        written = 0;
        for frameIndex = 1:frameCount
            if cancelExport || isClosing || ~isgraphics(fig), break; end
            sample = frameSamples(frameIndex);
            Render_Frame();
            drawnow;
            if cancelExport || isClosing || ~isgraphics(fig), break; end
            exportgraphics(layout, frameFile, 'Units', 'pixels', 'Width', 1000, ...
                'Height', imageHeight, 'Padding', 8, 'BackgroundColor', 'white');
            [indexed, palette] = rgb2ind(imread(frameFile), 256, 'nodither');
            if frameIndex == 1
                imwrite(indexed, palette, temporaryGIF, 'gif', ...
                    'LoopCount', Inf, 'DelayTime', delays(frameIndex));
            else
                imwrite(indexed, palette, temporaryGIF, 'gif', ...
                    'WriteMode', 'append', 'DelayTime', delays(frameIndex));
            end
            written = frameIndex;
            saveButton.String = sprintf('Cancel (%d%%)', round(100*frameIndex/frameCount));
        end
        if written == frameCount && ~cancelExport && ~isClosing
            [success, message] = movefile(temporaryGIF, file, 'f');
            assert(success, '%s', message);
            setappdata(fig, 'ReplayGIFFile', file);
            fprintf('Saved GIF: %s\n', file);
        end
        clear cleanup;
    end

    function Finish_Export(frameFile, temporaryGIF, savedSample, wasPlaying, savedResize, controls)
        if isfile(frameFile), delete(frameFile); end
        if isfile(temporaryGIF), delete(temporaryGIF); end
        isExporting = false;
        if ~isClosing && isgraphics(fig)
            fig.Resize = savedResize;
            set(controls, 'Enable', 'on');
            saveButton.String = 'Save GIF';
            Seek_Sample(savedSample);
            if wasPlaying && ~isClosing && isgraphics(fig)
                playButton.Value = 1;
                Toggle_Play();
            end
        end
    end

    function Close_Replay(varargin)
        isClosing = true;
        cancelExport = true;
        if isvalid(playbackTimer)
            stop(playbackTimer);
            delete(playbackTimer);
        end
        delete(fig);
    end
end

function value = Read_Base(name, varargin)
    if evalin('base', sprintf('exist(''%s'', ''var'')', name))
        value = evalin('base', name);
    elseif ~isempty(varargin)
        value = varargin{1};
    else
        error('BeamformingReplay:MissingVariable', ...
            'Missing workspace variable "%s". Run the dynamic script and this method first.', name);
    end
end
