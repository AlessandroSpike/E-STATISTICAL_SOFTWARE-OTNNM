function A = strength_driven_attachment(n, m, alpha)
    % Strength-Driven Attachment Model
    % This function generates a weighted undirected network using a strength-driven preferential attachment mechanism
    % Inputs:
    %   n: final number of nodes
    %   m: number of edges to attach from a new node to existing nodes
    %   alpha: preferential attachment exponent
    % Output:
    %   A: adjacency matrix of the generated network
    
    % Initialize the network with m+1 nodes
    A = zeros(n);
    % Create initial random connections between the first m+1 nodes
    A(1:m+1, 1:m+1) = triu(rand(m+1), 1);
    % Make the network undirected by adding the transpose
    A = A + A';
    
    % Add remaining nodes one by one
    for i = m+2:n
        % Calculate node strengths (sum of edge weights) for existing nodes
        strengths = sum(A(1:i-1, 1:i-1));
        
        % Calculate attachment probabilities based on node strengths
        probs = strengths.^alpha;
        probs = probs / sum(probs);  % Normalize to create probability distribution
        
        % Select target nodes for new connections
        targets = zeros(1, m);
        for j = 1:m
            % Use roulette wheel selection to choose a target node
            r = rand();
            cumprobs = cumsum(probs);
            targets(j) = find(cumprobs > r, 1);
            
            % Prevent duplicate connections by setting probability to 0
            probs(targets(j)) = 0;
            probs = probs / sum(probs);  % Renormalize probabilities
        end
        
        % Add new edges with random weights
        for j = 1:m
            weight = rand();  % Generate a random weight between 0 and 1
            A(i, targets(j)) = weight;
            A(targets(j), i) = weight;  % Ensure the network remains undirected
        end
    end
end

% Strength-Driven Attachment Model:
% This model extends the Barabási-Albert (BA) model to weighted networks, 
% considering node strengths (sum of edge weights) instead of degrees for 
% preferential attachment (Barrat et al., 2004).
% This code generates networks where the probability of attachment is 
% proportional to node strengths raised to power alpha. 
% This allows for tuning the preferential attachment strength and can 
% lead to different strength distributions.