%PLOT_RLS_RESULTS Plot RLS error, weights, constellation, and radiation pattern.
%   y, e: output/error; W: weight history; constellation: reference symbols.
%   phi: azimuth [rad]; patternRLS: normalized field.
%   phiDS, phiIS: source azimuths [rad]; postStart: first constellation sample.
%   mseWindow: averaging length; wo: Wiener weights on the same scale as W.
%   patternFloor_dB: radial minimum (default: -35 dB); SER: optional symbol error rate.

function figRLS = Plot_RLS_Results(y, e, W, constellation, phi, patternRLS, ...
    phiDS, phiIS, postStart, mseWindow, wo, patternFloor_dB, SER)
    if nargin < 12
        patternFloor_dB = -35;
    end
    if nargin < 13
        SER = NaN;
    end
    figRLS = Plot_LMS_Results(y, e, W, constellation, phi, patternRLS, ...
        phiDS, phiIS, postStart, mseWindow, wo, patternFloor_dB, SER, 'RLS');
end
