% Strength-Driven Attachment Network Study
%
% This script conducts a comprehensive study on strength-driven attachment networks.
% It explores how various parameters affect network properties and compares the original
% network with null models generated using Optimal Transport (OT).
%
% The study follows these main steps:
% 1) Generate synthetic networks using a strength-driven attachment model
% 2) Compute network statistics (in/out degree and strength)
% 3) Apply Optimal Transport to create a transport plan
% 4) Generate null models based on the OT plan
% 5) Compare properties of the original network with the null models
%
% The script uses parallel processing to speed up computations across multiple parameter combinations.
% Results are saved for further analysis, allowing for insights into how the strength-driven model
% and OT-based null models differ in their network properties.

clc; close all; clearvars
rng(1) % Set random seed for reproducibility

%% Parameters
% Each parameter is defined as a vector to explore different values
N = [100 500];        % Number of nodes (small and medium-sized networks)
M = [3 5];            % Number of links drawn for each new node
Gamma = [.01 .1];     % Entropy coefficient for Optimal Transport
Ripetition = [50 200]; % Number of null models to generate

ParamNames = {"Number of Nodes"; "Link Drawn"; "Entropy Parameter"; "Number of Riperitions"};
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
    m = Parame(2);        % Number of links drawn
    gamma = Parame(3);    % Entropy parameter for OT
    ripetition = Parame(4); % Number of null models to generate

    % Generate strength-driven attachment network
    A = strength_driven_attachment(n, m, 1);
    A = A / sum(sum(A));  % Normalize adjacency matrix

    % Compute network properties
    [in_strength, out_strength, in_degree, out_degree, rec_degree] = compute_network_properties(A);

    % Store original network statistics
    RealIndeg{x} = in_degree;
    RealOutdeg{x} = out_degree;
    RealInstr{x} = in_strength;
    RealOutstr{x} = out_strength;

    % Prepare for Optimal Transport
    instrprob = in_strength / sum(in_strength);
    outstrprob = out_strength / sum(out_strength);

    % Create cost matrix for OT
    C = out_degree * (in_degree');
    C = 1 ./ C;  % Inverse cost matrix

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
        NullW = sampled_joint;
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
save("strength_driven_attachment_results.mat", "PermOutstr", "PermInstr", ...
    "PermOutdeg", "PermIndeg", "RealOutstr", "RealInstr", "RealOutdeg", "RealIndeg", ...
    "ParamCombinations", "ParamNames")