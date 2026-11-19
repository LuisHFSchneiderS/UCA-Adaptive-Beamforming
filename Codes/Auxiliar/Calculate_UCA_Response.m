%CALCULATE_UCA_RESPONSE Compute the coupled receive response at a fixed azimuth.
%   K: element count; phiTX: arrival azimuth in the array plane [rad].
%   lambda, L, b, R: wavelength, dipole length, wire radius, and array radius [m].
%   ZL: K-by-K diagonal matrix of identical receiver loads [ohm].
%   h: K-by-1 complex response vector at the loaded antenna ports.

function h = Calculate_UCA_Response(K, lambda, L, b, R, phiTX, ZL)
    beta = 2*pi/lambda;

    S = Calculate_UCA_Elements_Distance_Matrix(K, R);
    Cpos = Define_UCA_Elements_Position(K, R);
    Z = Calculate_UCA_Impedance_Matrix(L, b, lambda, S);

    ZA = Z(1, 1);
    ZT = ZL(1, 1);

    steering_vector = zeros(K, 1);
    for i = 1:K
        xi = Cpos(i, 1);
        yi = Cpos(i, 2);

        steering_vector(i) = ...
            (lambda/(pi*sin(pi*L/lambda))) * (1-cos(pi*L/lambda)) * ...
            exp(-1i*beta*(xi*cos(phiTX) + yi*sin(phiTX)));
    end

    Cpl = (Z + ZL) \ ((ZT + ZA)*eye(K));

    h = (ZT/(ZT + ZA))*Cpl*steering_vector;
end