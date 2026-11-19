%ANIMATE_WIENER_RESULTS Replay workspace outputs with fixed Wiener weights.
%   Animate_Wiener_Results('Duration', 20, 'WindowSize', 2000)
%   Options: Duration [s], FPS, WindowSize [symbols], AutoPlay, GIFFile.

function fig = Animate_Wiener_Results(varargin)
    fig = Animate_Beamforming_Results('Wiener', varargin{:});
end
