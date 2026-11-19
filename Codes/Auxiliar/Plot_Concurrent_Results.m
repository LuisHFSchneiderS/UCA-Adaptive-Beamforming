%PLOT_CONCURRENT_RESULTS Plot concurrent CMA-DD beamforming results.
%   y: output; eMSE: offline IQ error; W: weight history.
%   constellation: reference symbols; phi: azimuth [rad].
%   patternConcurrent: normalized field; phiDS, phiIS: source azimuths [rad].
%   postStart: first constellation sample; mseWindow: averaging length.
%   patternFloor_dB: radial display minimum (default: -35 dB).
%   SER: optional final symbol error rate after offline correction.

function figConcurrent = Plot_Concurrent_Results(y, eMSE, W, constellation, phi, patternConcurrent, ...
    phiDS, phiIS, postStart, mseWindow, patternFloor_dB, SER)
    if nargin < 11
        patternFloor_dB = -35;
    end
    if nargin < 12
        SER = NaN;
    end
    figConcurrent = Plot_CMA_Results(y, eMSE, W, constellation, phi, patternConcurrent, ...
        phiDS, phiIS, postStart, mseWindow, 'Concurrent', patternFloor_dB, SER);
end
