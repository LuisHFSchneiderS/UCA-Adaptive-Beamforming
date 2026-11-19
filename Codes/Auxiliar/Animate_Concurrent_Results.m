%ANIMATE_CONCURRENT_RESULTS Replay workspace Concurrent results after phase correction.
%   Animate_Concurrent_Results('Duration', 20, 'WindowSize', 2000)
%   Options: Duration [s], FPS, WindowSize [symbols], AutoPlay, GIFFile.

function fig = Animate_Concurrent_Results(varargin)
    fig = Animate_Beamforming_Results('Concurrent', varargin{:});
end
