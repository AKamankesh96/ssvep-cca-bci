function [A, B, r, U, V] = cca_ssvep(X, Y)

% CCA_SSVEP Perform Canonical Correlation Analysis between EEG and reference signals.
%
% Inputs:
%   X - [nSamples x nChannels] EEG segment
%   Y - [nSamples x nComponents] reference signal matrix
%
% Outputs:
%   A, B - canonical coefficients
%   r - canonical correlations
%   U, V - canonical variates

    [n, p1] = size(X);
    p2 = size(Y, 2);

    if size(Y, 1) ~= n
        error('X and Y must have the same number of samples.');
    end
    if n <= 1
        error('At least two samples are required.');
    end

    % Mean-center both inputs
    X = X - mean(X, 1);
    Y = Y - mean(Y, 1);

    % QR factorization with rank handling for X
    [Q1, T11, perm1] = qr(X, 0);
    rankX = sum(abs(diag(T11)) > eps(abs(T11(1))) * max(n, p1));
    if rankX == 0
        error('Input X has rank 0.');
    elseif rankX < p1
        Q1  = Q1(:, 1:rankX);
        T11 = T11(1:rankX, 1:rankX);
    end

    % QR factorization with rank handling for Y
    [Q2, T22, perm2] = qr(Y, 0);
    rankY = sum(abs(diag(T22)) > eps(abs(T22(1))) * max(n, p2));
    if rankY == 0
        error('Input Y has rank 0.');
    elseif rankY < p2
        Q2  = Q2(:, 1:rankY);
        T22 = T22(1:rankY, 1:rankY);
    end

    % Compute canonical coefficients and correlations
    d = min(rankX, rankY);
    [L, D, M] = svd(Q1' * Q2, 0);

    A = T11 \ L(:, 1:d) * sqrt(n - 1);
    B = T22 \ M(:, 1:d) * sqrt(n - 1);
    r = min(max(diag(D(:, 1:d))', 0), 1);

    % Restore original variable order and dimensions
    A(perm1, :) = [A; zeros(p1 - rankX, d)];
    B(perm2, :) = [B; zeros(p2 - rankY, d)];

    % Return canonical variates if requested
    if nargout > 3
        U = X * A;
        V = Y * B;
    end
end