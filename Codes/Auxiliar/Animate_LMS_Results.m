%ANIMATE_LMS_RESULTS Replay workspace LMS results without running the filter.
%   Animate_LMS_Results('Duration', 20, 'WindowSize', 2000)
%   Options: Duration [s], FPS, WindowSize [symbols], AutoPlay, GIFFile.

function fig = Animate_LMS_Results(varargin)
    fig = Animate_Beamforming_Results('LMS', varargin{:});
end
