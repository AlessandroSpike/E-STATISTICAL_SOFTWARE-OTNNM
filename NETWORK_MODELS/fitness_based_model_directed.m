function A = fitness_based_model_directed(n, p, fitness_func, weight_func)
    % This function generates a directed weighted network based on node fitness
    % Inputs:
    %   n: number of nodes
    %   p: base probability of edge formation
    %   fitness_func: function handle for generating node fitness
    %   weight_func: function handle for generating edge weights
    % Output:
    %   A: adjacency matrix of the directed network

    % Initialize the adjacency matrix
    A = zeros(n);

    % Generate fitness values for all nodes
    fitness = fitness_func(n);
    
    % Iterate through all possible edges
    for i = 1:n
        for j = 1:n
            if i ~= j  % Avoid self-loops
                % Determine if an edge should be formed based on node fitness
                if rand < p * fitness(i) * fitness(j)
                    % Generate a weight for the edge
                    weight = weight_func();
                    % Assign the weight to the adjacency matrix
                    A(i, j) = weight;
                    % Note: A(j, i) is not set, making the network directed
                end
            end
        end
    end
end

% Example usage:
% A = fitness_based_model_directed(100, 0.1, @(n) rand(1, n), @() exprnd(1));

% Directed Fitness-Based Model:
% This model creates a directed network where the probability of edge
% formation from node i to node j depends on the fitness of both nodes.
% The resulting adjacency matrix A is not symmetric, representing
% directed edges. A(i,j) represents an edge from node i to node j.