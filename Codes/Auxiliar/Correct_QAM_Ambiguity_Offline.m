%CORRECT_QAM_AMBIGUITY_OFFLINE Select a fixed quadrant from the frame prefix.
%   symbols, reference: detected and transmitted symbol indices.
%   M: QAM order; evaluation: logical mask or first SER evaluation symbol.
%   yRaw: optional IQ samples, rotated before demodulation.
%   calibrationSymbols: number of initial reference symbols (default: up to 1000).
%   symbolsCorrected: corrected indices; rotation_deg: 0, 90, 180, or 270.
%   SER_before, SER_after: evaluation-window symbol error rates.
%   yCorrected: corrected IQ samples, when yRaw is supplied.

function [symbolsCorrected, rotation_deg, SER_before, SER_after, yCorrected] = ...
    Correct_QAM_Ambiguity_Offline(symbols, reference, M, evaluation, yRaw, calibrationSymbols)
    symbols = symbols(:);
    reference = reference(:);
    if nargin < 4
        evaluation = 1;
    end
    assert(sqrt(M) == fix(sqrt(M)), 'Quadrant correction requires square QAM.');
    assert(numel(symbols) == numel(reference), 'Symbol sequence lengths must match.');
    if islogical(evaluation)
        assert(isvector(evaluation) && numel(evaluation) == numel(symbols), ...
            'The SER evaluation mask must have one entry per symbol.');
        evaluationIndex = evaluation(:);
    else
        validateattributes(evaluation, {'numeric'}, ...
            {'scalar', 'integer', '>=', 1, '<=', numel(symbols)});
        evaluationIndex = evaluation:numel(symbols);
    end
    if nargin < 6
        calibrationSymbols = min(1000, numel(symbols));
    end
    validateattributes(calibrationSymbols, {'numeric'}, ...
        {'scalar', 'integer', '>=', 1, '<=', numel(symbols)});

    calibrationIndex = 1:calibrationSymbols;
    hasIQ = nargin >= 5;
    yCorrected = [];
    if hasIQ
        yRaw = yRaw(:);
        assert(numel(yRaw) == numel(symbols), 'IQ and symbol sequence lengths must match.');
        symbols = qamdemod(yRaw, M, 'UnitAveragePower', true);
    else
        alphabet = qammod((0:M-1).', M, 'UnitAveragePower', true);
    end
    SER_before = mean(symbols(evaluationIndex) ~= reference(evaluationIndex));

    symbolMaps = zeros(M, 4);
    calibrationSER = zeros(4, 1);

    for k = 1:4
        if hasIQ
            candidate = qamdemod((1i)^(k-1)*yRaw(calibrationIndex), M, ...
                'UnitAveragePower', true);
        else
            symbolMaps(:, k) = qamdemod((1i)^(k-1)*alphabet, M, ...
                'UnitAveragePower', true);
            candidate = symbolMaps(symbols(calibrationIndex)+1, k);
        end
        calibrationSER(k) = mean(candidate ~= reference(calibrationIndex));
    end

    [~, bestRotation] = min(calibrationSER);
    rotation_deg = 90*(bestRotation-1);
    if hasIQ
        yCorrected = (1i)^(bestRotation-1)*yRaw;
        symbolsCorrected = qamdemod(yCorrected, M, 'UnitAveragePower', true);
    else
        symbolsCorrected = symbolMaps(symbols+1, bestRotation);
    end
    SER_after = mean(symbolsCorrected(evaluationIndex) ~= reference(evaluationIndex));
end
