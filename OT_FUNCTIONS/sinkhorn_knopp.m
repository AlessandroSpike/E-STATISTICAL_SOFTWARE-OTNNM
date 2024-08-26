function [plan, totalCost] = sinkhorn_knopp(marginal1, marginal2, costMatrix, gamma, maxIter, tol)
    % Sinkhorn-Knopp Algorithm for Optimal Transport
    % This function computes the optimal transport plan with entropy regularization using the Sinkhorn-Knopp algorithm
    % 
    % Inputs:
    %   marginal1: vector of non-negative numbers summing to 1 (first marginal)
    %   marginal2: vector of non-negative numbers summing to 1 (second marginal)
    %   costMatrix: matrix where costMatrix(i,j) represents the cost of moving
    %               one unit of mass from location i in marginal1 to location j in marginal2
    %   gamma: regularization parameter
    %   maxIter: maximum number of iterations
    %   tol: tolerance for convergence
    % 
    % Outputs:
    %   plan: matrix where plan(i,j) represents the amount of mass moved from
    %         location i in marginal1 to location j in marginal2
    %   totalCost: the total cost of the transport plan

    % Ensure the marginals are valid probability distributions
    assert(abs(sum(marginal1) - 1) < 1e-9, 'First marginal does not sum to 1');
    assert(abs(sum(marginal2) - 1) < 1e-9, 'Second marginal does not sum to 1');
    assert(length(marginal1) == length(marginal2), 'Marginals must have the same length');
    assert(all(size(costMatrix) == [length(marginal1), length(marginal2)]), 'Cost matrix dimensions do not match marginals');
    
    % Get the dimension of the problem
    n = length(marginal1);

    % Initialize variables
    K = exp(-costMatrix / gamma);  % Compute the Gibbs kernel
    u = ones(n, 1);  % Initialize u vector
    v = ones(n, 1);  % Initialize v vector
    
    % Sinkhorn-Knopp iteration
    for iter = 1:maxIter
        u_prev = u;
        v_prev = v;
        
        % Update u and v
        u = marginal1 ./ (K * v);
        v = marginal2 ./ (K' * u);
        
        % Check for convergence
        if max(abs(u - u_prev)) < tol && max(abs(v - v_prev)) < tol
            break;
        end
    end
    
    % Compute the optimal transport plan
    plan = diag(u) * K * diag(v);
    
    % Compute the total cost of the transport plan
    totalCost = sum(sum(plan .* costMatrix));
end