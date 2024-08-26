% Fitness-Based Directed Network Model Study
%
% This script conducts a comprehensive study on fitness-based directed network models.
% It explores how various parameters affect network properties and compares the original
% network with null models generated using Optimal Transport (OT).
%
% The study follows these main steps:
% 1) Generate synthetic networks using a fitness-based directed model
% 2) Compute network statistics (in/out degree and strength)
% 3) Apply Optimal Transport to create a transport plan
% 4) Generate null models based on the OT plan
% 5) Compare properties of the original network with the null models
%
% The script uses parallel processing to speed up computations across multiple parameter combinations.
% Results are saved for further analysis, allowing for insights into how the fitness-based model
% and OT-based null models differ in their network properties.

clc; close all; clearvars
rng(1) % Set random seed for reproducibility

%% Parameters
% Each parameter is defined as a vector to explore different values

% N: Number of nodes in the network
%    Explores both small and medium-sized networks
N = [100 500];

% M: Probability of attachment in the fitness-based model
%    Higher values create denser networks
M = [.2 .5];

% Gamma: Entropy regularization parameter for Optimal Transport
%        Lower values make OT more similar to exact transport
Gamma = [.01 .1];

% Ripetition: Number of null models to generate for each parameter combination
%             Higher values provide more robust statistical comparisons
Ripetition = [50 200];

% Names of parameters for result labeling
ParamNames = {"Number of Nodes"; "Link Probability"; "Entropy Parameter"; "Number of Riperitions"};

% Generate all combinations of parameters
ParamCombinations = table2array(combinations(N, M, Gamma, Ripetition));

%% Initialize containers for results
% These cells will store results for each parameter combination
RealIndeg = cell(size(ParamCombinations, 1), 1);
RealOutdeg = cell(size(ParamCombinations, 1), 1);
RealInstr = cell(size(ParamCombinations, 1), 1);
RealOutstr = cell(size(ParamCombinations, 1), 1);

PermIndeg = cell(size(ParamCombinations, 1), 1);
PermOutdeg = cell(size(ParamCombinations, 1), 1);
PermInstr = cell(size(ParamCombinations, 1), 1);
PermOutstr = cell(size(ParamCombinations, 1), 1);

%% Main parallel loop over parameter combinations
parfor x = 1:size(ParamCombinations, 1)
    % Extract parameters for this combination
    Parame = ParamCombinations(x, :);
    n = Parame(1);
    m = Parame(2);
    gamma = Parame(3);
    ripetition = Parame(4);

    % Generate "real network" using fitness-based directed model
    A = fitness_based_model_directed(n, m, @(n) rand(1, n), @() exprnd(1));
    A = A / sum(sum(A));  % Normalize adjacency matrix

    % Compute network statistics
    [in_strength, out_strength, in_degree, out_degree, rec_degree] = compute_network_properties(A);

    % Store original network properties
    RealIndeg{x} = in_degree;
    RealOutdeg{x} = out_degree;
    RealInstr{x} = in_strength;
    RealOutstr{x} = out_strength;

    % Prepare for Optimal Transport
    instrprob = in_strength / sum(in_strength);
    outstrprob = out_strength / sum(out_strength);

    % Compute cost matrix for OT
    C = out_degree * (in_degree');
    C = 1 ./ C;
    maximum = max(C(~isinf(C)));
    C(isinf(C)) = maximum;

    % Apply Optimal Transport
    [T, totalCost] = sinkhorn_knopp(outstrprob, instrprob, C, gamma, 1000, 10^-200);

    % Generate null models
    indegE = zeros(n, ripetition);
    outdegE = zeros(n, ripetition);
    instrE = zeros(n, ripetition);
    outstrE = zeros(n, ripetition);

    for t = 1:ripetition
        % Sample from the OT joint distribution
        [X, Y, indices, marginal_X, marginal_Y, sampled_joint] = sample_2d_distribution(T, nnz(A));
        NullW = sampled_joint;
        NullB = NullW;
        NullB(NullB > 0) = 1;

        % Compute statistics for null model
        indegE(:, t) = sum(NullB);
        outdegE(:, t) = sum(NullB, 2);
        instrE(:, t) = sum(NullW);
        outstrE(:, t) = sum(NullW, 2);
    end

    % Store null model properties
    PermIndeg{x} = indegE;
    PermOutdeg{x} = outdegE;
    PermInstr{x} = instrE;
    PermOutstr{x} = outstrE;
end

%% Save results
save("fitness_based_model_directed_results.mat", "PermOutstr", "PermInstr", ...
    "PermOutdeg", "PermIndeg", "RealOutstr", "RealInstr", "RealOutdeg", "RealIndeg", "ParamCombinations", "ParamNames")