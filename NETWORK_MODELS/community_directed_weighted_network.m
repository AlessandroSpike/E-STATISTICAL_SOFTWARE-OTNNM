function [W, communities] = community_directed_weighted_network(n, c, p_in, p_out, sparsity, asymmetry)
    % This function generates a directed weighted network with community structure
    % Inputs:
    %   n: number of nodes
    %   c: number of communities
    %   p_in: probability of intra-community edges
    %   p_out: probability of inter-community edges
    %   sparsity: desired sparsity of the network
    %   asymmetry: level of asymmetry in edge weights
    % Outputs:
    %   W: adjacency matrix of the network
    %   communities: vector of community assignments for each node

    % Assign nodes to communities more evenly
    communities = repelem(1:c, ceil(n/c));
    communities = communities(randperm(length(communities)));
    communities = communities(1:n);

    % Initialize adjacency matrix
    W = zeros(n);

    % Create edges based on p_in and p_out probabilities
    for i = 1:n
        for j = 1:n
            if i ~= j  % Avoid self-loops
                if communities(i) == communities(j)
                    if rand < p_in
                        W(i,j) = 1;
                    end
                else
                    if rand < p_out
                        W(i,j) = 1;
                    end
                end
            end
        end
    end

    % Adjust network to achieve desired sparsity
    num_edges = round(n * n * sparsity);
    edges = find(W > 0);
    if length(edges) > num_edges
        % Remove excess edges randomly
        remove = randsample(edges, length(edges) - num_edges);
        W(remove) = 0;
    elseif length(edges) < num_edges
        % Add missing edges randomly, avoiding self-loops
        zero_edges = find(W == 0);
        diagonal_indices = 1:n+1:n^2;
        diagonal_indices(diagonal_indices > numel(W)) = [];
        zero_edges = setdiff(zero_edges, diagonal_indices);
        add = randsample(zero_edges, num_edges - length(edges));
        W(add) = 1;
    end

    % Assign weights to edges
    for i = 1:n
        for j = 1:n
            if W(i,j) > 0
                 if communities(i) == communities(j)
                    W(i,j) = 0.1 + 0.9*rand;  % Intra-community weights: [0.1, 1]
                else
                    W(i,j) = 1 + 4*rand;  % Inter-community weights: [1, 5]
                end
            end
        end
    end

    % Introduce noise in weights (currently set to 0)
    noise_level = 0;
    W = W .* (1 + noise_level * (rand(size(W)) - 0.5));

    % Introduce asymmetric flow/connectivity patterns
    for i = 1:n
        for j = 1:n
            if W(i,j) > 0 && rand < asymmetry
                W(j,i) = W(i,j) * (1 - asymmetry + asymmetry * rand);
            end
        end
    end

    % Add a few strong connections that don't align with community structure
    num_strong = round(n * 0.1);  % 10% of nodes get strong misaligned connections
    for i = 1:num_strong
        source = randi(n);
        target = randi(n);
        while source == target || communities(source) == communities(target)
            target = randi(n);
        end
        W(source, target) = 2 + rand;  % Strong inter-community connection
    end

    % Normalize weights to maintain heterogeneous strength distributions
    for i = 1:n
        out_strength = sum(W(i,:));
        if out_strength > 0
            W(i,:) = W(i,:) / out_strength * n;
        end
    end
end