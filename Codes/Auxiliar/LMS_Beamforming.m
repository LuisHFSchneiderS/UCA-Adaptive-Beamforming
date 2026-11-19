%LMS_BEAMFORMING Apply supervised complex LMS beamforming.
%   V_in: N-by-K samples; DS: N-by-1 reference; step: LMS step.
%   train: initial training fraction (default: 1); blocksize: subframe length.
%   y: output; e: reference error; W: K-by-N weights before updates.
%   validationMask: symbols outside each training prefix.

function [y, e, W, validationMask] = LMS_Beamforming(V_in, DS, step, train, blocksize)
    [Dstream_size, K] = size(V_in);
    if nargin < 4
        train = 1;
    end
    if nargin < 5
        blocksize = Dstream_size;
    end
    w = complex(zeros(K, 1));
    W = complex(zeros(K, Dstream_size));
    e = zeros(Dstream_size, 1);
    y = zeros(Dstream_size, 1);
    trainsize = round(blocksize*train);
    validationMask = mod((0:Dstream_size-1).', blocksize) >= trainsize;

    for i = 1:Dstream_size
        x = V_in(i, :).';

        y(i) = w'*x;
        
        e(i) = (DS(i) - y(i));

        W(:, i) = w;

        if ~validationMask(i)
            w = w + step.*x.*conj(e(i));
        end
    end
end
