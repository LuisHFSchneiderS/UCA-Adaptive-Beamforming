%DEFINE_UCA_ELEMENTS_POSITION Define the Cartesian centers of a receiving UCA.
%
%   C = Define_UCA_Elements_Position(K, R)
%
%   Places K antenna element centers uniformly on a circle of radius R in
%   the xy-plane, with the UCA center at the coordinate origin.
%
%   Inputs
%       K - Number of array elements; integer greater than or equal to 2.
%       R - UCA radius [m]; positive scalar.
%
%   Output
%       C - K-by-3 real Cartesian position matrix [m].
%           Row i contains [x_i, y_i, z_i], with z_i = 0 for every element.
%           C contains coordinates, not mutual-coupling coefficients.
%
%   Geometry and indexing
%       phi_i = 2*pi*(i-1)/K, for i = 1,...,K.
%       C(i,:) = [R*cos(phi_i), R*sin(phi_i), 0].
%   Azimuth increases counterclockwise from +x toward +y. The first element
%   lies on +x; the circular endpoint is not duplicated.
%
%   The element ordering matches Calculate_UCA_Elements_Distance_Matrix
%   and Calculate_UCA_Response. Dipole lengths and orientations are defined
%   by the electromagnetic model, rather than by this position matrix.
%
%   See also Calculate_UCA_Elements_Distance_Matrix,
%            Calculate_UCA_Response, Calculate_Radiation_Pattern.

function C = Define_UCA_Elements_Position(K, R)
    C = zeros(K, 3);

    for i = 1:K
        phi = 2*pi*(i-1)/K;
        C(i, 1) = R*cos(phi);
        C(i, 2) = R*sin(phi);
    end
end