function prepare_offline_data(inputMatFile, outputDir)

% PREPARE_OFFLINE_DATA Extract offline SSVEP epochs and labels from raw EEG.
%
% Steps:
%   1. Load raw recording
%   2. Detect flicker and target marker segments
%   3. Extract 5-second EEG epochs
%   4. Bandpass filter each epoch
%   5. Remove per-channel DC offset
%   6. Save processed data and labels

    if nargin < 2
        error('Usage: prepare_offline_data(inputMatFile, outputDir)');
    end

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    S = load(inputMatFile);
    if ~isfield(S, 'data')
        error('Input file must contain a variable named "data".');
    end

    rawData = squeeze(S.data);

    % Configuration parameters (adjust as needed)
    markerChannel = 26;
    eegChannels = [6:14, 18:21];
    trimStartSamples = 42200;
    samplesPerEpoch = 5*512;   % 5 s at 512 Hz
    fs = 512;
    bandpassRangeHz = [4, 45];
    nClasses = 5;
    trialsPerClass = 20;

    % Remove initial non-task samples
    rawData = rawData(:, trimStartSamples:end);
    markerSignal = rawData(markerChannel, :);

    % Detect flicker-active (marker value 4) and target-transitio period (marker value 3)
    flickerIdx = find(markerSignal == 4);
    targetIdx  = find(markerSignal == 3);

    if isempty(flickerIdx)
        error('No samples with marker value 4 were found.');
    end
    if isempty(targetIdx)
        error('No samples with marker value 3 were found.');
    end

    [~, flickerEndIdx] = findContiguousSegments(flickerIdx);
    [targetStartIdx, ~] = findContiguousSegments(targetIdx);

    if any(targetStartIdx <= 1)
        error('Cannot read target labels from the sample preceding the target event.');
    end

    % Original label mapping: label = preceding marker value - 6
    labels = markerSignal(targetStartIdx - 1) - 6;
    labels = labels(:);

    if any(labels < 1 | labels > nClasses)
        error('Detected labels are outside the expected range 1 to %d.', nClasses);
    end

    if numel(flickerEndIdx) ~= numel(labels)
        error('Mismatch between number of flicker segments (%d) and labels (%d).', ...
            numel(flickerEndIdx), numel(labels));
    end

    % Epoch start index follows the original implementation
    epochStartIdx = flickerEndIdx - (samplesPerEpoch - 1);

    % Group epoch starts by class
    targetEpochStarts = cell(nClasses, 1);
    for classId = 1:nClasses
        targetEpochStarts{classId} = epochStartIdx(labels == classId);
    end

    for classId = 1:nClasses
        if numel(targetEpochStarts{classId}) < trialsPerClass
            error('Class %d has only %d trials; expected at least %d.', ...
                classId, numel(targetEpochStarts{classId}), trialsPerClass);
        end
    end

    % Extract epochs as [class x channel x sample x trial]
    eegByClass = zeros(nClasses, numel(eegChannels), samplesPerEpoch, trialsPerClass);

    for classId = 1:nClasses
        for trialId = 1:trialsPerClass
            s0 = targetEpochStarts{classId}(trialId);
            s1 = s0 + samplesPerEpoch - 1;

            if s0 < 1 || s1 > size(rawData, 2)
                error('Epoch exceeds data bounds for class %d, trial %d.', classId, trialId);
            end

            eegByClass(classId, :, :, trialId) = rawData(eegChannels, s0:s1);
        end
    end

    % Convert to [nTrials x nChannels x nSamples]
    data = [];
    for classId = 1:nClasses
        classData = permute(squeeze(eegByClass(classId, :, :, :)), [3, 1, 2]);
        data = cat(1, data, classData);
    end

    % Filter and mean-center each epoch
    [b, a] = butter(4, bandpassRangeHz / (fs / 2), 'bandpass');

    for trialId = 1:size(data, 1)
        epoch = squeeze(data(trialId, :, :));   % [channels x samples]
        epoch = filtfilt(b, a, epoch')';
        epoch = epoch - mean(epoch, 2);
        data(trialId, :, :) = epoch;
    end

    save(fullfile(outputDir, 'Data_Offline5.mat'), 'data', '-v7');
    save(fullfile(outputDir, 'Label_Offline5.mat'), 'labels', '-v7');

    fprintf('Saved processed data to:\n  %s\n  %s\n', ...
        fullfile(outputDir, 'Data_Offline5.mat'), ...
        fullfile(outputDir, 'Label_Offline5.mat'));
end

function [segmentStarts, segmentEnds] = findContiguousSegments(indexVector)
%FINDCONTIGUOUSSEGMENTS Return starts and ends of contiguous index runs.
    breaks = [true, diff(indexVector) > 1];
    segmentStarts = indexVector(breaks);
    segmentEnds = indexVector([find(breaks(2:end)) - 1, numel(indexVector)]);
end