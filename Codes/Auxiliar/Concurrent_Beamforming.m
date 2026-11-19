%CONCURRENT_BEAMFORMING Apply concurrent CMA and decision-directed beamforming.
%   V_in: N-by-K samples; stepCMA, stepDD: adaptive steps.
%   gamma: dispersion constant; constellation: reference symbols.
%   y: output; e: modulus error; W: weights before updates.

function [y, e, W] = Concurrent_Beamforming(V_in, stepCMA, stepDD, gamma, constellation)
    [Dstream_size, K] = size(V_in);
    constellation = constellation(:);

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

        wCMA = w + stepCMA*x*conj(y(i))*e(i);
        yCMA = wCMA'*x;

        [~, decision] = min(abs(constellation - y(i)));
        [~, decisionCMA] = min(abs(constellation - yCMA));

        w = wCMA;
        if decision == decisionCMA  % Consistent decisions enable DD.
            eDD = constellation(decision) - y(i);
            w = w + stepDD*x*conj(eDD);
        end
    end
end
