%PLOT_INTERFERENCE_RAYS Add translucent radial guides for interferers.
%   polarAxes: polar axes; phiIS: interferer azimuths [rad].
%   rays: handles to the radial guides.

function rays = Plot_Interference_Rays(polarAxes, phiIS)
    rays = gobjects(numel(phiIS), 1);
    if isempty(phiIS)
        return;
    end
    wasHeld = ishold(polarAxes);
    hold(polarAxes, 'on');
    for i = 1:numel(phiIS)
        rays(i) = polarplot(polarAxes, [phiIS(i) phiIS(i)], polarAxes.RLim, ...
            'Color', [0.45 0.45 0.45], 'LineWidth', 1.2, ...
            'HandleVisibility', 'off', 'HitTest', 'off', 'PickableParts', 'none', ...
            'Tag', 'InterferenceRay');
        rays(i).Color = [0.45 0.45 0.45 0.30];
    end
    uistack(rays, 'bottom');
    limitListener = addlistener(polarAxes, 'RLim', 'PostSet', ...
        @(~, ~) Update_Interference_Rays(polarAxes, rays));
    setappdata(polarAxes, 'InterferenceRayListener', limitListener);
    if ~wasHeld
        hold(polarAxes, 'off');
    end
end

function Update_Interference_Rays(polarAxes, rays)
    if ~isgraphics(polarAxes)
        return;
    end
    for i = 1:numel(rays)
        if isgraphics(rays(i))
            rays(i).RData = polarAxes.RLim;
        end
    end
end
