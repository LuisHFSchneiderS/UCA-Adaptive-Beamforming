%CALCULATE_UCA_IMPEDANCE_MATRIX Compute analytical UCA port impedances.
%
%   Z = Calculate_UCA_Impedance_Matrix(L, b, lambda, S)
%
%   Implements the Schelkunoff self- and mutual-impedance expressions used
%   in Professor Fernando de Castro's lecture notes. The model describes
%   identical, parallel, center-fed thin dipoles in free space, with their
%   centers in a plane perpendicular to the dipole axes. In this project,
%   the dipoles are vertical and their centers form a UCA in the xy-plane.
%
%   Inputs
%       L      - Total physical length of each receiving dipole [m].
%       b      - Wire radius of each receiving dipole [m].
%       lambda - Carrier wavelength [m].
%       S      - K-by-K center-to-center distance matrix [m].
%                S must be real and symmetric, with zero diagonal and
%                positive off-diagonal entries. K is inferred from size(S,1).
%
%   Output
%       Z      - K-by-K complex antenna port impedance matrix [ohm].
%                Z(i,i) is the self-impedance of element i.
%                Z(i,j), i ~= j, is the mutual impedance between i and j.
%                Reciprocity gives Z(i,j) = Z(j,i), without conjugation.
%
%   Model conventions and scope
%       - All dipoles share the same length L and wire radius b.
%       - L is the full length; the internal variable h equals L/2.
%       - L need not equal lambda/2, within the validity and nonsingular
%         parameter range of the thin-wire analytical approximation.
%       - The phasor time convention is exp(1j*omega*t).
%       - Free-space impedance is approximated as eta0 = 120*pi ohm.
%       - Receiver loads and external-wave excitation are not included.
%
%   For a receiver load matrix ZL and open-circuit excitation vector u,
%   the corresponding load-voltage vector is
%       v = ZL * ((Z + ZL) \ u).
%
%   Dependencies
%       Symbolic Math Toolbox: eulergamma, sinint, and cosint.
%
%   Reference
%       F. C. C. de Castro, Comunicacoes Estrategicas, Chapters II.1-II.5,
%       Schelkunoff impedance routines, slides 144-148.
%
%   See also Calculate_UCA_Elements_Distance_Matrix, Calculate_UCA_Response.

function Z = Calculate_UCA_Impedance_Matrix(L, b, lambda, S)
    [K, ~] = size(S);
    em = double(eulergamma);
    
    beta = 2*pi/lambda;
    h = L/2;
    q = beta*h;
    Zc = 120*(log(L/b)-1);
    F = em + log(2*q) - cosint(2*q);

    Rstar = 60*(F + 0.5*sin(2*q)*(sinint(4*q)-2*sinint(2*q)) + 0.5*cos(2*q)*(em+log(q)+cosint(4*q)-2*cosint(2*q)));
    Xstar = 60*(sinint(2*q) - 0.5*sin(2*q)*(em+log(q)-cosint(4*q)) - (0.5*cos(2*q)*sinint(4*q)));

    M = 60*(F - 1 + cos(2*q));
    N = 60*(sinint(2*q) - sin(2*q));
    W = (1 + M/Zc)*sin(q) + (Xstar + N -1i*Rstar)*cos(q)/Zc;

    Z = zeros(K);

    for i = 1:K
        for j = 1:K
            if i == j % self-impedance
                Z(i,j) = (Rstar*sin(q) + 1i*((Xstar - N)*sin(q) - (Zc - M)*cos(q)))/W;
            else % mutual impedance
                d = S(i, j);

                A1 = beta*d;
                A2 = beta*sqrt(d^2 + h^2) + q;
                A3 = A2 -2*q;
                A4 = beta*sqrt(d^2 + L^2) - 2*q;
                A5 = A4 + 4*q;

                S1 = sinint(A1);
                S2 = sinint(A2);
                S3 = sinint(A3);
                S4 = sinint(A4);
                S5 = sinint(A5);

                C1 = cosint(A1);
                C2 = cosint(A2);
                C3 = cosint(A3);
                C4 = cosint(A4);
                C5 = cosint(A5);

                O1 = 2*C1 - 2*(C2 + C3) + C4 + C5;
                O2 = 4*C1 - 2*(C2 + C3);
                O3 = 2*(S3 - S2) + S5 - S4;

                Q1 = 2*(S2 + S3) - 2*S1 - S4 - S5;
                Q2 = 2*(S2 + S3) - 4*S1;
                Q3 = 2*(C3 - C2) - C4 +C5;

                BR = O1*cos(2*q)+O2+O3*sin(2*q);
                BX = Q1*cos(2*q)+Q2+Q3*sin(2*q);

                Z(i, j) = ((30)/(W^2))*(BR + 1i*BX);
            end
        end
    end
end