%PLOT_CMA_SDD_RESULTS Plot CMA-SDD MSE, weights, constellation, and pattern.
%   Same inputs as Plot_Concurrent_Results. SER is an optional input.

function figCMASDD = Plot_CMA_SDD_Results(y, eMSE, W, constellation, phi, patternCMASDD, ...
    phiDS, phiIS, postStart, mseWindow, patternFloor_dB, SER)
    if nargin < 11
        patternFloor_dB = -35;
    end
    if nargin < 12
        SER = NaN;
    end
    figCMASDD = Plot_CMA_Results(y, eMSE, W, constellation, phi, patternCMASDD, ...
        phiDS, phiIS, postStart, mseWindow, 'CMA-SDD', patternFloor_dB, SER);
end
