%CALCULATE_RADIATION_PATTERN Evaluate a complex far-field azimuth cut.
%
%   E = Calculate_Radiation_Pattern(I0, r, phi_numel, L, C, lambda)
%
%   Evaluates the theta-polarized electric-field phasor for identical,
%   parallel, z-directed thin dipoles at a fixed observation distance r.
%   The observation zenith angle is theta = pi/2, so the returned samples
%   describe the horizontal (H-plane) cut over a full azimuth revolution.
%
%   Inputs
%       I0        - Vector of K complex radiation-current amplitudes [A],
%                   using peak phasors. Each value is the coefficient of
%                   the corresponding dipole's sinusoidal current profile.
%                   For L = lambda/2, it equals the feed-current phasor.
%       r         - Observation distance from the coordinate origin [m];
%                   positive and in the far field of the complete array.
%       phi_numel - Number of azimuth samples; integer greater than or
%                   equal to 2. Both 0 and 2*pi are included.
%       L         - Common total physical dipole length [m]; positive scalar.
%       C         - K-by-3 Cartesian position matrix [m]. Row k is the
%                   center [x_k, y_k, z_k] of dipole k. K = size(C,1).
%                   This C denotes positions, not a coupling matrix.
%       lambda    - Carrier wavelength [m]; positive scalar.
%
%   Output
%       E         - phi_numel-by-1 complex peak E_theta field samples [V/m],
%                   ordered according to phi = linspace(0,2*pi,phi_numel).'.
%                   The output is a field phasor, without normalization
%                   or conversion to power, gain, or decibels.
%
%   Intended field model at theta = pi/2
%       beta = 2*pi/lambda;
%       E(phi) = (60/r)*exp(1j*(pi/2-beta*r))*(1-cos(pi*L/lambda))
%                * sum_k I0(k)*exp(1j*beta*(x_k*cos(phi)+y_k*sin(phi))).
%   The element factor equals one for half-wave dipoles. The z-coordinate
%   contributes zero propagation phase in this horizontal far-field cut.
%   The phasor time convention is exp(1j*omega*t), and free-space impedance
%   is approximated as eta0 = 120*pi ohm.
%
%   Coupling and excitation
%   Supply radiation currents that already account for mutual coupling,
%   obtained from the coupled impedance model or an equivalent transmit
%   excitation for the receive beamformer. Digital combining weights require
%   conversion to equivalent currents before being supplied as I0.
%
%   For a nonzero field, a normalized azimuth pattern can be formed as
%       F = abs(E)/max(abs(E));
%       FdB = 20*log10(F);
%   Use the same azimuth grid as the field samples when plotting.
%
%   Reference
%       F. C. C. de Castro, Comunicacoes Estrategicas, Chapters II.1-II.5,
%       slide 68, equation (14), specialized to theta = pi/2.
%
%   See also Define_UCA_Elements_Position, Calculate_UCA_Impedance_Matrix.

function E = Calculate_Radiation_Pattern(I0, r, phi_numel, L, C, lambda)
    [K, ~] = size(C);
    beta = 2*pi/lambda;

    phi = linspace(0, 2*pi, phi_numel).';

    E = zeros(phi_numel, 1);
    for i = 1:K
        Cx = C(i, 1);
        Cy = C(i, 2);
        Ei = I0(i).*exp(1i*beta*(Cx*cos(phi) + Cy*sin(phi))).*(1 - cos(pi*L/lambda))*(60/r)*exp(1i*((pi/2)-beta*r));
        E = E+Ei;
    end
end