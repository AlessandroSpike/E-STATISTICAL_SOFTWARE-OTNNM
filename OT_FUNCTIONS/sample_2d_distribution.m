function [X, Y, indices, marginal_X, marginal_Y, sampled_joint] = sample_2d_distribution(P, k)
    % Sample from a 2D Probability Distribution
    % This function samples points from a given 2D probability distribution
    % Inputs:
    %   P: A 2-D matrix representing the joint probability distribution
    %   k: Number of samples to draw (default is 1 if not provided)
    % Outputs:
    %   X: Vector of sampled x-coordinates
    %   Y: Vector of sampled y-coordinates
    %   indices: Vector of linear indices corresponding to the sampled points
    %   marginal_X: Empirical marginal distribution along X-axis
    %   marginal_Y: Empirical marginal distribution along Y-axis
    %   sampled_joint: Empirical joint distribution from samples

    % Set default value for k if not provided
    if nargin < 2
        k = 1;
    end

    % Ensure P is a valid probability distribution
    if any(P(:) < 0) || abs(sum(P(:)) - 1) > 1e-10
        error('Input matrix P must be a valid probability distribution');
    end

    % Flatten the 2-D matrix into a 1-D array
    P_flat = P(:);

    % Compute the cumulative sum for inverse transform sampling
    cumsum_P = cumsum(P_flat);

    % Initialize output vectors
    X = zeros(1, k);
    Y = zeros(1, k);
    indices = zeros(1, k);

    % Sample k points using inverse transform sampling
    for i = 1:k
        % Generate a random number
        r = rand();

        % Find the index where the random number falls
        indices(i) = find(cumsum_P >= r, 1, 'first');

        % Convert the linear index back to 2-D coordinates
        [Y(i), X(i)] = ind2sub(size(P), indices(i));
    end

    % Compute empirical marginal distributions from samples
    [n, m] = size(P);
    marginal_X = histcounts(X, 1:m+1) / k;
    marginal_Y = histcounts(Y, 1:n+1) / k;
    
    % Ensure marginal_Y is a column vector to match the original distribution's shape
    marginal_Y = marginal_Y(:);

    % Compute sampled 2-D joint distribution
    sampled_joint = zeros(size(P));
    for i = 1:k
        sampled_joint(Y(i), X(i)) = sampled_joint(Y(i), X(i)) + 1;
    end
    sampled_joint = sampled_joint / k;
end