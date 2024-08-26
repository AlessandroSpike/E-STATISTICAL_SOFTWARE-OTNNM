% Weighted Random Network Study
%
% This script conducts a comprehensive study on weighted random directed networks.
% It explores how various parameters affect network properties and compares the original
% network with null models generated using Optimal Transport (OT).
%
% The study follows these main steps:
% 1) Generate synthetic weighted directed networks
% 2) Compute network statistics (in/out degree and strength)
% 3) Apply Optimal Transport to create a transport plan
% 4) Generate null models based on the OT plan
% 5) Compare properties of the original network with the null models
%
% The script uses parallel processing to speed up computations across multiple parameter combinations.
% Results are saved for further analysis, allowing for insights into how the original model
% and OT-based null models differ in their network properties.

clc; close all; clearvars
rng(1) % Set random seed for reproducibility

%% Parameters
% Each parameter is defined as a vector to explore different values
N = [100 500];        % Number of nodes (small and medium-sized networks)
M = [.1 .3];          % Probability of link creation (network density)
Gamma = [.01 .1];     % Entropy coefficient for Optimal Transport
Ripetition = [50 200]; % Number of null models to generate

% Create descriptive names for parameters and generate all combinations
ParamNames = {"Number of Nodes"; "Link Probability"; "Entropy Parameter"; "Number of Riperitions"};
ParamCombinations = table2array(combinations(N, M, Gamma, Ripetition));

%% Initialize result containers
% These cells will store results for each parameter combination
% 'Real' prefix for original network, 'Perm' prefix for null models
RealIndeg = cell(size(ParamCombinations,1), 1);
RealOutdeg = cell(size(ParamCombinations,1), 1);
RealInstr = cell(size(ParamCombinations,1), 1);
RealOutstr = cell(size(ParamCombinations,1), 1);

PermIndeg = cell(size(ParamCombinations,1), 1);
PermOutdeg = cell(size(ParamCombinations,1), 1);
PermInstr = cell(size(ParamCombinations,1), 1);
PermOutstr = cell(size(ParamCombinations,1), 1);

%% Main parallel loop
% This loop processes each parameter combination in parallel
parfor x = 1:size(ParamCombinations, 1)
    Parame = ParamCombinations(x,:);
    n = Parame(1);        % Number of nodes
    m = Parame(2);        % Link probability
    gamma = Parame(3);    % Entropy parameter for OT
    ripetition = Parame(4); % Number of null models to generate

    % Generate weighted directed graph
    A = weighted_directed_graph(n, m, @() rand(1));
    A = A / sum(sum(A));  % Normalize adjacency matrix to create a probability matrix

    % Compute network properties of the original network
    [in_strength, out_strength, in_degree, out_degree, rec_degree] = compute_network_properties(A);

    % Store original network statistics
    RealIndeg{x} = in_degree;
    RealOutdeg{x} = out_degree;
    RealInstr{x} = in_strength;
    RealOutstr{x} = out_strength;

    % Prepare for Optimal Transport
    instrprob = in_strength / sum(in_strength);  % Normalize in-strength to probability
    outstrprob = out_strength / sum(out_strength);  % Normalize out-strength to probability
    C = out_degree * (in_degree');  % Create cost matrix
    C = 1 ./ C;  % Inverse cost matrix for OT

    % Perform Optimal Transport using Sinkhorn-Knopp algorithm
    [T, totalCost] = sinkhorn_knopp(outstrprob, instrprob, C, gamma, 1000, 1e-200);

    % Generate null models
    indegE = zeros(n, ripetition);
    outdegE = zeros(n, ripetition);
    instrE = zeros(n, ripetition);
    outstrE = zeros(n, ripetition);

    for t = 1:ripetition
        % Sample from joint distribution obtained from OT
        [X, Y, indices, marginal_X, marginal_Y, sampled_joint] = sample_2d_distribution(T, nnz(A));
        NullW = sampled_joint;  % Weighted null model
        NullB = NullW > 0;  % Binary version of null model

        % Compute null model statistics
        indegE(:,t) = sum(NullB);  % In-degree
        outdegE(:,t) = sum(NullB, 2);  % Out-degree
        instrE(:,t) = sum(NullW);  % In-strength
        outstrE(:,t) = sum(NullW, 2);  % Out-strength
    end

    % Store null model statistics
    PermIndeg{x} = indegE;
    PermOutdeg{x} = outdegE;
    PermInstr{x} = instrE;
    PermOutstr{x} = outstrE;
end

%% Save results
% Store all computed statistics and parameters for further analysis
save("weighted_random_graph_results.mat", "PermOutstr", "PermInstr", ...
    "PermOutdeg", "PermIndeg", "RealOutstr", "RealInstr", "RealOutdeg", "RealIndeg", ...
    "ParamCombinations", "ParamNames")