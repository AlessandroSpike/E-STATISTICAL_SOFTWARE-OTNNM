% Community Detection in Directed Weighted Networks Study
%
% This script conducts a comprehensive study on community detection in directed weighted networks.
% It explores the performance of two community detection methods:
% 1) Standard modularity optimization
% 2) Optimal Transport (OT) based modularity optimization
%
% The study generates synthetic networks with known community structure using various parameters.
% It then applies both detection methods and compares their performance using two metrics:
% a) Variation of Information (VI) between detected and true communities
% b) Modularity score of the detected community structure
%
% The script uses parallel processing to speed up computations across multiple parameter combinations.
% Results are saved for further analysis, allowing for insights into how network properties and
% detection methods influence community detection accuracy in directed weighted networks.

clc; close all; clearvars
rng(1) % Set random seed for reproducibility

%% Parameters
% Each parameter is defined as a vector to explore different values

% K: Number of communities in the network
%    Higher values create more complex community structures
K = [5 10]; 

% P_IN: Probability of intra-community connections
%       Higher values create more densely connected communities
P_IN = [.6 .8]; 

% P_OUT: Probability of inter-community connections
%        Higher values create more connections between communities
P_OUT = [.2 .4]; 

% Gamma: Entropy regularization parameter for Optimal Transport
%        Lower values make OT more similar to exact transport
Gamma = [.00001 .0001]; 

% Sparsity: Overall network sparsity
%           Higher values create denser networks
Sparsity = [.1 .8];

% Asymmetry: Degree of asymmetry in edge weights
%            Higher values create more asymmetric connections
Asymmetry = [.1 .8];

% Number of repetitions for each parameter combination
Riperition = 50;

% Names of parameters for result labeling
ParamNames = {"Number of Communities"; "Within Link Prob."; "Between Link Prob."; "Entropy Parameter"; ...
    "Sparsity"; "Asymmetry"};

% Generate all combinations of parameters
ParamCombinations = table2array(combinations(K, P_IN, P_OUT, Gamma, Sparsity, Asymmetry));

%% Initialize containers for results
% These arrays will store results for each parameter combination and repetition
RealVariation = zeros(size(ParamCombinations, 1), Riperition);
PermVariation = zeros(size(ParamCombinations, 1), Riperition);
RealModul = zeros(size(ParamCombinations, 1), Riperition);
PermModul = zeros(size(ParamCombinations, 1), Riperition);

%% Main parallel loop over parameter combinations
for x = 1:size(ParamCombinations, 1)
    % Extract parameters for this combination
    Parame = ParamCombinations(x, :);
    k = Parame(1);
    p_in = Parame(2);
    p_out = Parame(3);
    gamma = Parame(4);
    sparsity = Parame(5);
    asymmetry = Parame(6);
    ripetition = Riperition;

    parfor t = 1:ripetition
        % Generate "real network" with known community structure
        [A, communities] = community_directed_weighted_network(200, k, p_in, p_out, sparsity, asymmetry);
        A = A / sum(sum(A));  % Normalize adjacency matrix

        % Compute network statistics
        [in_strength, out_strength, in_degree, out_degree, rec_degree] = compute_network_properties(A);
        instrprob = in_strength / sum(in_strength);
        outstrprob = out_strength / sum(out_strength);

        % Optimal transport part
        C = out_degree * (in_degree');  % Cost matrix
        C = 1 ./ C;  % Inverse cost
        [T, totalCost] = sinkhorn_knopp(outstrprob, instrprob, C, gamma, 1000, 10^-200);

        % Find communities using standard modularity
        [Ci, Q] = modularity_dir(A, 1);
        % Find communities using OT-based modularity
        [CiOT, QOT] = OTmodularity_dir(A, 1, T);

        % Store modularity results
        RealModul(x, t) = Q;
        PermModul(x, t) = QOT;

        % Compute variation of information for standard modularity
        VIn = variation_of_information(Ci, communities');
        RealVariation(x, t) = VIn;

        % Compute variation of information for OT-based modularity
        VInOT = variation_of_information(CiOT, communities');
        PermVariation(x, t) = VInOT;
    end
end

%% Save results
save("community_results.mat", "PermVariation", ...
    "RealVariation", "RealModul", "PermModul", "ParamCombinations", "ParamNames")