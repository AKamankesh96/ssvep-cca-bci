function references = generate_ssvep_references(frequencies, phases, fs, durationSec, nHarmonics)

% GENERATE_SSVEP_REFERENCES Generate sine/cosine reference templates for SSVEP.
%
% Inputs:
%   frequencies - target frequencies in Hz
%   phases - phase offsets in radians
%   fs - sampling frequency in Hz
%   durationSec - window duration in seconds
%   nHarmonics  - number of harmonics to include
%
% Output:
%   references - cell array of reference matrices [2*nHarmonics x nSamples]

    if nargin < 5
        nHarmonics = 3;
    end

    if numel(frequencies) ~= numel(phases)
        error('frequencies and phases must have the same number of elements.');
    end

    nTargets = numel(frequencies);
    nSamples = round(durationSec * fs);
    t = (0:nSamples - 1) / fs;

    references = cell(nTargets, 1);

    for k = 1:nTargets
        ref = zeros(2 * nHarmonics, nSamples);
        row = 1;

        for h = 1:nHarmonics
            ref(row, :)     = sin(2 * pi * h * frequencies(k) * t + phases(k));
            ref(row + 1, :) = cos(2 * pi * h * frequencies(k) * t + phases(k));
            row = row + 2;
        end

        references{k} = ref;
    end
end