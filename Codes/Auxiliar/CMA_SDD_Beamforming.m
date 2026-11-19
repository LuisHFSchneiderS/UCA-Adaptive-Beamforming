%CMA_SDD_BEAMFORMING Apply concurrent CMA and soft decision-directed beamforming.
%   V_in: N-by-K samples; stepCMA, stepSDD: adaptive steps.
%   gamma: dispersion constant; constellation: square QAM symbols.
%   rho: per-axis cluster variance, in squared constellation units.
%   y, e, W: output, modulus error, and weights before updates.

function [y, e, W] = CMA_SDD_Beamforming(V_in, stepCMA, stepSDD, gamma, constellation, rho)
    [Dstream_size, K] = size(V_in);
    constellation = constellation(:);
    realLevels = unique(real(constellation));
    imagLevels = unique(imag(constellation));
    realCenters = (realLevels(1:2:end) + realLevels(2:2:end))/2;
    imagCenters = (imagLevels(1:2:end) + imagLevels(2:2:end))/2;

    w = complex(zeros(K, 1));
    w(1) = 1;
    W = complex(zeros(K, Dstream_size));
    e = zeros(Dstream_size, 1);
    y = complex(zeros(Dstream_size, 1));

    for i = 1:Dstream_size
        x = V_in(i, :).';
        y(i) = w'*x;
        e(i) = gamma - abs(y(i))^2;
        W(:, i) = w;

        [~, rBlock] = min(abs(realCenters - real(y(i))));
        [~, iBlock] = min(abs(imagCenters - imag(y(i))));
        rPair = realLevels(2*rBlock-1:2*rBlock);
        iPair = imagLevels(2*iBlock-1:2*iBlock);
        symbols = [rPair(1) + 1i*iPair(1); rPair(1) + 1i*iPair(2); ...
            rPair(2) + 1i*iPair(1); rPair(2) + 1i*iPair(2)];
        distances = abs(y(i) - symbols).^2;
        probabilities = exp(-(distances-min(distances))/(2*rho));
        probabilities = probabilities/sum(probabilities);
        softDecision = sum(probabilities.*symbols);
        eSDD = softDecision - y(i);

        w = w + stepCMA*x*conj(y(i))*e(i) + stepSDD*x*conj(eSDD);
    end
end
