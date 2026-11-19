%RLS_BEAMFORMING Apply supervised complex RLS beamforming.
%   V_in: N-by-K samples; DS: N-by-1 reference.
%   lambdaRLS: forgetting factor (0, 1]; delta: regularization, P(0) = I/delta.
%   train: initial training fraction (default: 1); blocksize: subframe length.
%   y: output; e: reference error; W: K-by-N weights before updates.
%   validationMask: symbols outside each training prefix.

function [y, e, W, validationMask] = RLS_Beamforming(V_in, DS, lambdaRLS, delta, train, blocksize)
    [Dstream_size, K] = size(V_in);
    if nargin < 5
        train = 1;
    end
    if nargin < 6
        blocksize = Dstream_size;
    end
    w = complex(zeros(K, 1));
    W = complex(zeros(K, Dstream_size));
    e = complex(zeros(Dstream_size, 1));
    y = zeros(Dstream_size, 1);
    P = complex(eye(K)/delta);
    trainsize = round(blocksize*train);
    validationMask = mod((0:Dstream_size-1).', blocksize) >= trainsize;

    for i = 1:Dstream_size
        x = V_in(i, :).';

        y(i) = w'*x;
        
        e(i) = (DS(i) - y(i));

        W(:, i) = w;

        if ~validationMask(i)
            q = x'*P*x;

            g = (P*x)/(lambdaRLS + q);

            w = w + g*conj(e(i));

            P = (P - g*x'*P)/lambdaRLS;
        end
    end
end
