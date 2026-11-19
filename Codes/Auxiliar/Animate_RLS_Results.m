%ANIMATE_RLS_RESULTS Replay workspace RLS results without running the filter.
%   Animate_RLS_Results('Duration', 20, 'WindowSize', 2000)
%   Options: Duration [s], FPS, WindowSize [symbols], AutoPlay, GIFFile.

function fig = Animate_RLS_Results(varargin)
    fig = Animate_Beamforming_Results('RLS', varargin{:});
end
