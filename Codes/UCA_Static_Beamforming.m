clc; clear; close all;

codeFolder = fileparts(mfilename('fullpath'));  % Script folder.
if ~isfolder(fullfile(codeFolder, 'Auxiliar'))
    codeFolder = fileparts(matlab.desktop.editor.getActiveFilename);
end
addpath(fullfile(codeFolder, 'Auxiliar'));

%% General Parameters
f = 850e6;                       % Frequency [Hz].
c = physconst('LightSpeed');     % Light speed [m/s].
lambda = c/f;                   % Wavelength [m].

%% RX UCA Design
K = 6;                   % Dipole count.
L = 0.5*lambda;          % Dipole length [m].
I0 = ones(K, 1);         % Radiation currents [A, peak].
b = 0.005;              % Dipole wire radius [m].
R = lambda/4;           % UCA radius [m].
r = 10000;              % Observation distance [m].
phi_numel = 721;         % Azimuth sample count.
patternFloor_dB = -30;   % Normalized polar plot minimum [dB].

C = Define_UCA_Elements_Position(K, R);                 % Positions [m].
S = Calculate_UCA_Elements_Distance_Matrix(K, R);        % Distances [m].
Z = Calculate_UCA_Impedance_Matrix(L, b, lambda, S);     % Impedances [ohm].
elementSpacing = 2*R*sin(pi/K);  % Adjacent spacing [m].

Ifeed = sin(pi*L/lambda)*I0(:);  % Feed currents [A, peak].

Pacc = 0.5*real(Ifeed'*Z*Ifeed);  % Accepted power [W].

phi = linspace(0, 2*pi, phi_numel).';  % Azimuth [rad].
E = Calculate_Radiation_Pattern(I0, r, phi_numel, L, C, lambda);  % Field [V/m].
Eiso = sqrt(60*Pacc)/r;  % Isotropic field [V/m].

G = (abs(E)/Eiso).^2;  % Linear gain.
G_dBi = 10*log10(G);   % Gain [dBi].

%% RX UCA Gain and Geometry
showStructure = 'on';       % Structure: 'on' or 'off'.
structureFraction = 0.3;   % Relative structure radius.

% figUCA = Plot_UCA_Gain(phi, G_dBi, C, f, lambda, L, b, R, Pacc, showStructure, structureFraction);

%% TX Desired Signal
phiDS = 0*pi/180;
Dstream_size = 1e6;  % Symbol count.
Dmod_order = 16;         % Modulation order.

Dstream = randi([0 Dmod_order-1], Dstream_size, 1);  % Symbol indices.

DS = qammod(Dstream, Dmod_order, 'UnitAveragePower', true);  % Normalized IQ.

%% TX Interferent Signals
IS_config = [
    % Azimuth [deg]   QAM order
          20              16;
         120              16;
         180              16;
         240              16;
];
if isempty(IS_config)
    IS_config = zeros(0, 2);
end
IS_config = array2table(IS_config, 'VariableNames', ...
    {'Azimuth_deg', 'QAM_Order'});

Num_IS = height(IS_config);                    % Interferer count.
phiIS = IS_config.Azimuth_deg*pi/180;           % Azimuth [rad].
Istream_size = Dstream_size;  % Symbol count.

Istream = zeros(Istream_size, Num_IS);           % Symbol indices.
IS = complex(zeros(Istream_size, Num_IS));       % Interferer IQ.
for i = 1:Num_IS
    Istream(:, i) = randi([0 IS_config.QAM_Order(i)-1], Istream_size, 1);
    IS(:, i) = qammod(Istream(:, i), IS_config.QAM_Order(i), ...
        'UnitAveragePower', true);
end

%% RX Noise and SINR
SNR_RX = 25;           % Desired-signal SNR [dB].
SINR_target_dB = 10;   % Input SINR target with interferers [dB].

%% UCA elements terminals
ZL = 50*eye(K);

hDS = Calculate_UCA_Response(K, lambda, L, b, R, phiDS, ZL);

hIS = zeros(K, Num_IS);
for i = 1:Num_IS
    hIS(:, i) = Calculate_UCA_Response(K, lambda, L, b, R, phiIS(i), ZL);
end

%% Interference power adjustment
powerDS = mean(abs(DS).^2)*mean(abs(hDS).^2);
noisePower = powerDS/10^(SNR_RX/10);

powerIS = 0;
gainIS = zeros(Num_IS, 1);  % IQ amplitude gains.
targetPowerIS = 0;
targetPowerPerIS = 0;
if Num_IS > 0
    targetPowerIS = powerDS/10^(SINR_target_dB/10) - noisePower;
    targetPowerPerIS = targetPowerIS/Num_IS;

    for i = 1:Num_IS
        currentPowerIS = mean(abs(IS(:, i)).^2)*mean(abs(hIS(:, i)).^2);
        gainIS(i) = sqrt(targetPowerPerIS/currentPowerIS);
        IS(:, i) = IS(:, i)*gainIS(i);
        powerIS = powerIS + currentPowerIS*gainIS(i)^2;
    end
end

SINR_in = powerDS/(powerIS + noisePower);
SINR_in_dB = 10*log10(SINR_in);  % Nominal input SINR [dB].

%% Received signals and noise
V_in = DS*hDS.';
for i = 1:Num_IS
    V_in = V_in + IS(:, i)*hIS(:, i).';
end

noiseRX = sqrt(noisePower/2)*( ...
    randn(size(V_in)) + 1i*randn(size(V_in)));

V_in = V_in + noiseRX;

%% Training and validation
train_until = 0.2;      % Training fraction.
blocksize = Dstream_size;
trainsize = round(blocksize*train_until);
trainingMask = mod((0:Dstream_size-1).', blocksize) < trainsize;
validationMask = ~trainingMask;

%% Wiener solution
inputScale = sqrt(mean(abs(V_in).^2, 'all'));
[yWiener, eWiener, woWiener, Rx_est, pxd_est, JminWiener] = ...
    Wiener_Beamforming(V_in/inputScale, DS, trainingMask);

%% Wiener Results
mseWindow = 1000;  % MSE averaging window.
DstreamWiener = qamdemod(yWiener, Dmod_order, 'UnitAveragePower', true);
SER_Wiener = mean(DstreamWiener(validationMask) ~= Dstream(validationMask));
wWiener = woWiener/inputScale;          % Load-voltage weights.

IfeedWiener = (Z + ZL)'\(ZL'*wWiener);  % Equivalent feed currents.
I0Wiener = IfeedWiener/sin(pi*L/lambda);
I0Wiener = I0Wiener/max(abs(I0Wiener));
E_Wiener = Calculate_Radiation_Pattern(I0Wiener, r, phi_numel, L, C, lambda);
patternWiener = abs(E_Wiener)/max(abs(E_Wiener));  % Normalized field pattern.

constellation = qammod((0:Dmod_order-1).', Dmod_order, ...
    'UnitAveragePower', true);          % Reference constellation.

figWiener = Plot_Wiener_Results(yWiener, eWiener, woWiener, JminWiener, constellation, ...
    phi, patternWiener, phiDS, phiIS(1:Num_IS), mseWindow, patternFloor_dB, SER_Wiener);

%% LMS solution
stepLMS = 0.01;
[yLMS, eLMS, WLMS, validationMaskLMS] = ...
    LMS_Beamforming(V_in/inputScale, DS, stepLMS, train_until, blocksize);

%% LMS Results
postStartLMS = floor(0.8*Dstream_size) + 1;  % Post-transient start.
DstreamLMS = qamdemod(yLMS, Dmod_order, 'UnitAveragePower', true);
SER_LMS = mean(DstreamLMS(validationMaskLMS) ~= Dstream(validationMaskLMS));
wLMS = WLMS(:, end)/inputScale;            % Last stored weights.

IfeedLMS = (Z + ZL)'\(ZL'*wLMS);  % Equivalent feed currents.
I0LMS = IfeedLMS/sin(pi*L/lambda);        % Radiation currents.
I0LMS = I0LMS/max(abs(I0LMS));            % Relative current amplitudes.
E_LMS = Calculate_Radiation_Pattern(I0LMS, r, phi_numel, L, C, lambda);
patternLMS = abs(E_LMS)/max(abs(E_LMS));  % Normalized field pattern.

figLMS = Plot_LMS_Results(yLMS, eLMS, WLMS, constellation, phi, patternLMS, ...
    phiDS, phiIS(1:Num_IS), postStartLMS, mseWindow, woWiener, patternFloor_dB, SER_LMS);

%% RLS solution
lambdaRLS = 1;     % Forgetting factor.
deltaRLS = 1e-2;   % Initial regularization.
[yRLS, eRLS, WRLS, validationMaskRLS] = ...
    RLS_Beamforming(V_in/inputScale, DS, lambdaRLS, deltaRLS, train_until, blocksize);

%% RLS Results
postStartRLS = floor(0.8*Dstream_size) + 1;  % Post-transient start.
DstreamRLS = qamdemod(yRLS, Dmod_order, 'UnitAveragePower', true);
SER_RLS = mean(DstreamRLS(validationMaskRLS) ~= Dstream(validationMaskRLS));
wRLS = WRLS(:, end)/inputScale;            % Last stored weights.

IfeedRLS = (Z + ZL)'\(ZL'*wRLS);  % Equivalent feed currents.
I0RLS = IfeedRLS/sin(pi*L/lambda);  % Radiation currents.
I0RLS = I0RLS/max(abs(I0RLS));      % Relative current amplitudes.
E_RLS = Calculate_Radiation_Pattern(I0RLS, r, phi_numel, L, C, lambda);
patternRLS = abs(E_RLS)/max(abs(E_RLS));  % Normalized field pattern.

figRLS = Plot_RLS_Results(yRLS, eRLS, WRLS, constellation, phi, patternRLS, ...
    phiDS, phiIS(1:Num_IS), postStartRLS, mseWindow, woWiener, patternFloor_dB, SER_RLS);

%% LMS and RLS comparison
figLMSRLSComparison = Plot_LMS_RLS_Comparison(yLMS, yRLS, eLMS, eRLS, constellation, ...
    phi, patternLMS, patternRLS, phiDS, phiIS(1:Num_IS), ...
    max(postStartLMS, postStartRLS), mseWindow, patternFloor_dB, [SER_LMS SER_RLS]);

%% CMA solution
stepCMA = 0.0002;
gammaCMA = mean(abs(constellation).^4)/mean(abs(constellation).^2);
[yCMA, eCMA, WCMA] = CMA_Beamforming(V_in/inputScale, stepCMA, gammaCMA);

%% CMA Results
postStartCMA = floor(0.8*Dstream_size) + 1;  % Post-transient start.
wCMA = WCMA(:, end)/inputScale;            % Last stored weights.

IfeedCMA = (Z + ZL)'\(ZL'*wCMA);    % Equivalent feed currents.
I0CMA = IfeedCMA/sin(pi*L/lambda);   % Radiation currents.
I0CMA = I0CMA/max(abs(I0CMA));       % Relative current amplitudes.
E_CMA = Calculate_Radiation_Pattern(I0CMA, r, phi_numel, L, C, lambda);
patternCMA = abs(E_CMA)/max(abs(E_CMA));  % Normalized field pattern.

%% Concurrent solution
stepConcurrentCMA = 0.0002;  % CMA step.
stepConcurrentDD = 0.0008;   % DD step.
gammaConcurrent = mean(abs(constellation).^4)/mean(abs(constellation).^2);
[yConcurrent, eConcurrent, WConcurrent] = Concurrent_Beamforming(V_in/inputScale, ...
    stepConcurrentCMA, stepConcurrentDD, gammaConcurrent, constellation);

%% Concurrent Results
postStartConcurrent = floor(0.8*Dstream_size) + 1;  % Post-transient start.
wConcurrent = WConcurrent(:, end)/inputScale;          % Last stored weights.

IfeedConcurrent = (Z + ZL)'\(ZL'*wConcurrent);  % Equivalent feed currents.
I0Concurrent = IfeedConcurrent/sin(pi*L/lambda); % Radiation currents.
I0Concurrent = I0Concurrent/max(abs(I0Concurrent));  % Relative current amplitudes.
E_Concurrent = Calculate_Radiation_Pattern(I0Concurrent, r, phi_numel, L, C, lambda);
patternConcurrent = abs(E_Concurrent)/max(abs(E_Concurrent));  % Normalized field pattern.

%% CMA-SDD solution
stepCMASDDCMA = 0.0002;  % CMA step.
stepCMASDDSDD = 0.0008;  % SDD step.
rhoCMASDD = 0.4/(2*(Dmod_order-1)/3);  % Unit-power cluster variance.
gammaCMASDD = mean(abs(constellation).^4)/mean(abs(constellation).^2);
[yCMASDD, eCMASDD, WCMASDD] = CMA_SDD_Beamforming(V_in/inputScale, ...
    stepCMASDDCMA, stepCMASDDSDD, gammaCMASDD, constellation, rhoCMASDD);

%% CMA-SDD Results
postStartCMASDD = floor(0.8*Dstream_size) + 1;  % Post-transient start.
wCMASDD = WCMASDD(:, end)/inputScale;  % Last stored weights.

IfeedCMASDD = (Z + ZL)'\(ZL'*wCMASDD);  % Equivalent feed currents.
I0CMASDD = IfeedCMASDD/sin(pi*L/lambda);
I0CMASDD = I0CMASDD/max(abs(I0CMASDD));
E_CMASDD = Calculate_Radiation_Pattern(I0CMASDD, r, phi_numel, L, C, lambda);
patternCMASDD = abs(E_CMASDD)/max(abs(E_CMASDD));  % Normalized field pattern.

%% Offline demodulation
DstreamCMA = qamdemod(yCMA, Dmod_order, 'UnitAveragePower', true);
DstreamConcurrent = qamdemod(yConcurrent, Dmod_order, 'UnitAveragePower', true);
DstreamCMASDD = qamdemod(yCMASDD, Dmod_order, 'UnitAveragePower', true);

%% Offline phase ambiguity correction
ambiguitySymbols = min(1000, trainsize);  % Initial reference symbols.
[DstreamCMA_corrected, rotationCMA_deg, SER_CMA_before, SER_CMA, yCMA_corrected] = ...
    Correct_QAM_Ambiguity_Offline(DstreamCMA, Dstream, Dmod_order, validationMask, yCMA, ambiguitySymbols);
[DstreamConcurrent_corrected, rotationConcurrent_deg, SER_Concurrent_before, SER_Concurrent, yConcurrent_corrected] = ...
    Correct_QAM_Ambiguity_Offline(DstreamConcurrent, Dstream, Dmod_order, validationMask, yConcurrent, ambiguitySymbols);
[DstreamCMASDD_corrected, rotationCMASDD_deg, SER_CMASDD_before, SER_CMASDD, yCMASDD_corrected] = ...
    Correct_QAM_Ambiguity_Offline(DstreamCMASDD, Dstream, Dmod_order, validationMask, yCMASDD, ambiguitySymbols);

%% Offline MSE
eMSE_CMA = DS - yCMA_corrected;       % IQ error after phase correction.
eMSE_Concurrent = DS - yConcurrent_corrected;
eMSE_CMASDD = DS - yCMASDD_corrected;
MSE_CMA = mean(abs(eMSE_CMA(postStartCMA:end)).^2);
MSE_Concurrent = mean(abs(eMSE_Concurrent(postStartConcurrent:end)).^2);
MSE_CMASDD = mean(abs(eMSE_CMASDD(postStartCMASDD:end)).^2);

%% CMA and Concurrent plots
figCMA = Plot_CMA_Results(yCMA_corrected, eMSE_CMA, WCMA, constellation, phi, ...
    patternCMA, phiDS, phiIS(1:Num_IS), postStartCMA, mseWindow, 'CMA', patternFloor_dB, SER_CMA);
figConcurrent = Plot_Concurrent_Results(yConcurrent_corrected, eMSE_Concurrent, WConcurrent, constellation, ...
    phi, patternConcurrent, phiDS, phiIS(1:Num_IS), postStartConcurrent, mseWindow, patternFloor_dB, SER_Concurrent);

%% CMA-SDD plot
figCMASDD = Plot_CMA_SDD_Results(yCMASDD_corrected, eMSE_CMASDD, WCMASDD, constellation, ...
    phi, patternCMASDD, phiDS, phiIS(1:Num_IS), postStartCMASDD, mseWindow, patternFloor_dB, SER_CMASDD);

%% CMA and Concurrent comparison
figComparison = Plot_CMA_Concurrent_Comparison(yCMA_corrected, yConcurrent_corrected, ...
    eMSE_CMA, eMSE_Concurrent, constellation, phi, patternCMA, patternConcurrent, ...
    phiDS, phiIS(1:Num_IS), max(postStartCMA, postStartConcurrent), mseWindow, ...
    patternFloor_dB, [SER_CMA SER_Concurrent]);

%% Concurrent and CMA-SDD comparison
figConcurrentCMASDDComparison = Plot_Concurrent_CMA_SDD_Comparison( ...
    yConcurrent_corrected, yCMASDD_corrected, eMSE_Concurrent, eMSE_CMASDD, ...
    constellation, phi, patternConcurrent, patternCMASDD, phiDS, phiIS(1:Num_IS), ...
    max(postStartConcurrent, postStartCMASDD), mseWindow, ...
    patternFloor_dB, [SER_Concurrent SER_CMASDD]);
