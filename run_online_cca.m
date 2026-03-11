function predictedLabels = run_online_cca(dataFile, labelFile)

% RUN_ONLINE_CCA Simulate online SSVEP target detection using fixed buffers.
%
% Inputs:
%   dataFile  - path to data file containing variable "data"
%   labelFile - optional path to label file containing variable "labels"
%
% Output:
%   predictedLabels - predicted class for each full 1-second buffer

    if nargin < 1
        error('Usage: run_online_cca(dataFile, labelFile)');
    end

    Sdata = load(dataFile);
    if ~isfield(Sdata, 'data')
        error('Data file must contain variable "data".');
    end
    data = Sdata.data;

    hasLabels = false;
    if nargin >= 2 && ~isempty(labelFile)
        Slabel = load(labelFile);
        if isfield(Slabel, 'labels')
            labels = Slabel.labels; %#ok<NASGU>
            hasLabels = true;
        end
    end

    % Configuration parameters (adjust as needed)
    fs = 512;
    windowSec = 1;
    windowSize = fs * windowSec;
    channels = 1:13;
    frequencies = [6.67, 7.5, 8.57, 10, 12];
    phases = [0, 0.5*pi, 0.5*pi, 0, pi];
    nHarmonics = 3;
    nTargets = numel(frequencies);

    % Generate reference templates
    references = generate_ssvep_references(frequencies, phases, fs, windowSec, nHarmonics);

    % Convert to [samples x channels x trials]
    data = permute(data, [3, 2, 1]);

    % Use the first trial for online simulation
    dataOnline = data(1:windowSize * 5, :, 1);

    buffer = zeros(0, size(dataOnline, 2));
    predictedLabels = [];

    for sampleIdx = 1:size(dataOnline, 1)
        buffer = [buffer; dataOnline(sampleIdx, :)]; %#ok<AGROW>

        if size(buffer, 1) == windowSize
            X = buffer(:, channels);
            X = X - mean(X, 1);

            rho = zeros(1, nTargets);
            for targetId = 1:nTargets
                [~, ~, r] = cca_ssvep(X, references{targetId}');
                rho(targetId) = max(r);
            end

            [~, predictedLabel] = max(rho);
            predictedLabels(end + 1, 1) = predictedLabel; %#ok<AGROW>
            fprintf('Detected target: %d\n', predictedLabel);

            % Reset for the next non-overlapping window
            buffer = zeros(0, size(dataOnline, 2));
        end
    end

    % Optional label note
    if hasLabels
        fprintf('Predictions generated for %d windows.\n', numel(predictedLabels));
        fprintf(['Labels were loaded, but evaluation depends on how the ', ...
                 'streamed windows are aligned with the ground truth.\n']);
    end
end