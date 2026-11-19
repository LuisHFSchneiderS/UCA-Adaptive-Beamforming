%CALCULATE_UCA_DYNAMIC_RESPONSE Compute coupled receive responses for a moving source.
%   lenDS: symbol count; K: element count; SR: symbol rate [symbols/s].
%   lambda, L, b, R: wavelength, dipole length, wire radius, and array radius [m].
%   phiTX_0: first-symbol azimuth [rad]; ang_vel: constant angular velocity [rad/s].
%   ZL: K-by-K diagonal matrix of identical receiver loads [ohm].
%   h: K-by-lenDS complex matrix of loaded-port responses, one column per symbol.

function h = Calculate_UCA_Dynamic_Response(lenDS, K, lambda, L, b, R, phiTX_0, ang_vel, SR, ZL)
    beta = 2*pi/lambda;
    
    n = 0:lenDS-1;

    phiTX = (ang_vel/SR)*n + phiTX_0;

    S = Calculate_UCA_Elements_Distance_Matrix(K, R);
    Cpos = Define_UCA_Elements_Position(K, R);
    Z = Calculate_UCA_Impedance_Matrix(L, b, lambda, S);

    ZA = Z(1, 1);
    ZT = ZL(1, 1);

    steering_vector = zeros(K, lenDS);
    for i = 1:K
        xi = Cpos(i, 1);
        yi = Cpos(i, 2);

        for j = 1:lenDS
            steering_vector(i, j) = ...
                (lambda/(pi*sin(pi*L/lambda))) * (1-cos(pi*L/lambda)) * ...
                exp(-1i*beta*(xi*cos(phiTX(j)) + yi*sin(phiTX(j))));
        end
    end

    Cpl = (Z + ZL) \ ((ZT + ZA)*eye(K));

    h = (ZT/(ZT + ZA))*Cpl*steering_vector;
end