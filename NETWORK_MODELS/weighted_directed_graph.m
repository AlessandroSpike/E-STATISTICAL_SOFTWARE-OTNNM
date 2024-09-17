function W = weighted_directed_graph(n, p, w_func)
    % Weighted Directed Graph Model
    % This function generates a weighted directed graph using a probabilistic approach
    % Inputs:
    %   n: number of nodes
    %   p: probability of edge formation
    %   w_func: function handle for generating weights
    % Output:
    %   W: adjacency matrix of the generated weighted directed graph

    % Initialize the adjacency matrix with zeros
    W = zeros(n);

    % Iterate through all possible edges
    for i = 1:n
        for j = 1:n
            % Avoid self-loops (edges from a node to itself)
            if i ~= j
                % Decide whether to create an edge based on probability p
                if rand < p
                    % If an edge is created, assign a weight using w_func
                    W(i, j) = w_func();
                    % Note: A(j, i) is not set, making the graph directed
                end
            end
        end
    end
end

% Example usage:
% W = weighted_directed_graph(100, 0.1, @() exprnd(1));

% Weighted Directed Graph Model:
% This model extends the Erdős-Rényi (ER) model to weighted directed networks. 
% Edges are formed between nodes with probability p, and their 
% weights are assigned using the function w_func. This allows 
% for flexibility in the weight distribution, such as exponential, 
% uniform, or other custom distributions.