%ANIMATE_CMA_RESULTS Replay workspace CMA results after offline phase correction.
%   Animate_CMA_Results('Duration', 20, 'WindowSize', 2000)
%   Options: Duration [s], FPS, WindowSize [symbols], AutoPlay, GIFFile.

function fig = Animate_CMA_Results(varargin)
    fig = Animate_Beamforming_Results('CMA', varargin{:});
end
