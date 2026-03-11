function results = run_offline_cca(dataFile, labelFile)

% RUN_OFFLINE_CCA Perform offline SSVEP target detection using CCA.
%
% Inputs:
%   dataFile - path to processed data file containing variable "data"
%   labelFile - path to label file containing variable "labels"
%
% Output:
%   results - struct with predictions, confusion matrix, and accuracy metrics

    if nargin < 2
        error('Usage: run_offline_cca(dataFile, labelFile)');
    end

    tic;

    Sdata = load(dataFile);
    Slabel = load(labelFile);

    if ~isfield(Sdata, 'data')
        error('Data file must contain variable "data".');
    end
    if ~isfield(Slabel, 'labels')
        error('Label file must contain variable "labels".');
    end

    data = Sdata.data;
    labels = Slabel.labels;

   % Configuration parameters (adjust as needed)
    fs = 512;
    windowSec = 5;
    frequencies = [6.67, 7.5, 8.57, 10, 12];
    phases = [0, 0.5*pi, 0.5*pi, 0, pi];
    channels = 1:13;
    nHarmonics = 3;
    nTargets = numel(frequencies);
    windowSamples = windowSec * fs;

    % Generate reference templates
    references = generate_ssvep_references(frequencies, phases, fs, windowSec, nHarmonics);

    % Convert to [samples x channels x trials]
    data = permute(data, [3, 2, 1]);
    data = data(1:windowSamples, :, :);

    nTrials = size(data, 3);
    predLabels = zeros(nTrials, 1);

    % CCA-based target detection
    for trialId = 1:nTrials
        X = data(:, :, trialId);
        X = X - mean(X, 1);

        rho = zeros(1, nTargets);
        for targetId = 1:nTargets
            [~, ~, r] = cca_ssvep(X(:, channels), references{targetId}');
            rho(targetId) = max(r);
        end

        [~, predLabels(trialId)] = max(rho);
    end

    % Evaluate predictions
    confMat = confusionmat(labels, predLabels);
    totalAccuracy = 100 * sum(diag(confMat)) / sum(confMat(:));

    classAccuracy = zeros(nTargets, 1);
    for targetId = 1:nTargets
        classAccuracy(targetId) = 100 * confMat(targetId, targetId) / sum(confMat(targetId, :));
    end

    elapsedTime = toc;

    fprintf('Confusion matrix:\n');
    disp(confMat);
    fprintf('Total accuracy: %.2f%%\n', totalAccuracy);
    for targetId = 1:nTargets
        fprintf('Accuracy Target %d: %.2f%%\n', targetId, classAccuracy(targetId));
    end
    fprintf('Elapsed time: %.4f s\n', elapsedTime);

    results = struct();
    results.predictedLabels = predLabels;
    results.confusionMatrix = confMat;
    results.totalAccuracy = totalAccuracy;
    results.classAccuracy = classAccuracy;
    results.elapsedTimeSec = elapsedTime;
    results.frequenciesHz = frequencies;
    results.phasesRad = phases;
end