%CALCULATE_UCA_ELEMENTS_DISTANCE_MATRIX Compute pairwise UCA center distances.
%
%   S = Calculate_UCA_Elements_Distance_Matrix(K, R)
%
%   Returns the center-to-center Euclidean distances between K receiving
%   elements uniformly distributed on a circle of radius R in the xy-plane.
%   The distance matrix supplies the separations used by the analytical
%   self- and mutual-impedance model.
%
%   Inputs
%       K - Number of array elements; integer greater than or equal to 2.
%       R - UCA radius [m]; positive scalar.
%
%   Output
%       S - K-by-K real, nonnegative, symmetric distance matrix [m].
%           S(i,j) is the distance between element centers i and j.
%           S(i,i) = 0, and distinct centers have positive separation.
%
%   Geometry and indexing
%       phi_i = 2*pi*(i-1)/K, for i = 1,...,K.
%       p_i = [R*cos(phi_i), R*sin(phi_i), 0].
%   Azimuth increases counterclockwise from +x toward +y. The first element
%   lies on +x, and the UCA center is the coordinate origin.
%
%   Chord-distance expression
%       S(i,j) = 2*R*abs(sin(pi*(i-j)/K)).
%   S describes geometry only; it contains no propagation phase, antenna
%   impedance, or receiver load information.
%
%   See also Calculate_UCA_Impedance_Matrix, Calculate_UCA_Response.

function S = Calculate_UCA_Elements_Distance_Matrix(K, R)
    S = zeros(K);
    
    for i = 1:K
        for j = 1:K
            S(i, j) = 2*R*abs(sin(pi*(i-j)/K));
        end
    end   
end