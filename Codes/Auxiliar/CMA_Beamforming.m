%CMA_BEAMFORMING Apply blind complex CMA beamforming.
%   V_in: N-by-K samples; step: CMA step; gamma: dispersion constant.
%   y: output; e: modulus error; W: K-by-N weights before updates.

function [y, e, W] = CMA_Beamforming(V_in, step, gamma)
    [Dstream_size, K] = size(V_in);

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

        w = w + step*x*conj(y(i))*e(i);
    end
end
