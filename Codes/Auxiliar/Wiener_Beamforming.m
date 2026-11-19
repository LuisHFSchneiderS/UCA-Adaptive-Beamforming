%WIENER_BEAMFORMING Estimate a fixed complex Wiener beamformer.
%   V_in: N-by-K samples; DS: N-by-1 reference; trainingMask: samples used for fitting.
%   y, e: output/error; wo: K-by-1 fixed weights.
%   Rx_est, pxd_est: training correlations; Jmin: minimum training MSE.

function [y, e, wo, Rx_est, pxd_est, Jmin] = Wiener_Beamforming(V_in, DS, trainingMask)
    [Dstream_size, K] = size(V_in);
    if nargin < 3
        trainingMask = true(Dstream_size, 1);
    end
    trainingIndex = find(trainingMask);
    Ntrain = numel(trainingIndex);
    Rx_est = complex(zeros(K));
    pxd_est = complex(zeros(K, 1));

    for i = trainingIndex(:).'
        x = V_in(i, :).';
        Rx_est = Rx_est + x*x';
        pxd_est = pxd_est + x*conj(DS(i));
    end
    Rx_est = Rx_est/Ntrain;
    pxd_est = pxd_est/Ntrain;

    wo = Rx_est \ pxd_est;
    Jmin = real(mean(abs(DS(trainingMask)).^2) - pxd_est'*wo);

    y = V_in*conj(wo);
    e = DS - y;
end
