%CONFIGURE_POLAR_DB_AXES Set native dB limits and physical data tips.
%   polarAxes, levelLine: axes and line with RData already in dB.
%   limits_dB: radial limits; unit: 'dB' or 'dBi'.

function Configure_Polar_dB_Axes(polarAxes, levelLine, limits_dB, unit)
    validateattributes(limits_dB, {'numeric'}, {'real', 'finite', 'vector', 'numel', 2});
    assert(limits_dB(2) > limits_dB(1), 'Radial limits must be increasing.');
    unit = validatestring(unit, {'dB', 'dBi'});
    rawLevels_dB = levelLine.RData;
    setappdata(levelLine, 'RawLevel_dB', rawLevels_dB);

    polarAxes.RLim = limits_dB;
    polarAxes.RTickLabelMode = 'auto';
    polarAxes.RAxis.TickLabelFormat = ['%g ' unit];
    polarAxes.Clipping = 'on';
    levelLine.Clipping = 'on';
    levelLine.DataTipTemplate.DataTipRows = [ ...
        dataTipTextRow('Azimuth [deg]', levelLine.ThetaData*180/pi); ...
        dataTipTextRow(['Level [' unit ']'], rawLevels_dB)];
    Update_Polar_dB_Display(polarAxes, levelLine, rawLevels_dB);
    limitListener = addlistener(polarAxes, 'RLim', 'PostSet', ...
        @(~, ~) Update_Polar_dB_Display(polarAxes, levelLine, rawLevels_dB));
    setappdata(levelLine, 'RadialLimitListener', limitListener);
end

function Update_Polar_dB_Display(polarAxes, levelLine, rawLevels_dB)
    if ~isgraphics(polarAxes) || ~isgraphics(levelLine)
        return;
    end
    displayedLevels_dB = rawLevels_dB;
    displayedLevels_dB(rawLevels_dB < polarAxes.RLim(1)) = polarAxes.RLim(1);
    displayedLevels_dB(rawLevels_dB > polarAxes.RLim(2)) = polarAxes.RLim(2);
    levelLine.RData = displayedLevels_dB;
    limits_dB = polarAxes.RLim;
    tickStep_dB = 5*max(1, ceil(diff(limits_dB)/50));
    polarAxes.RTick = unique([limits_dB(1), ...
        ceil(limits_dB(1)/tickStep_dB)*tickStep_dB:tickStep_dB:limits_dB(2), limits_dB(2)]);
    polarAxes.RTickLabelMode = 'auto';
end
