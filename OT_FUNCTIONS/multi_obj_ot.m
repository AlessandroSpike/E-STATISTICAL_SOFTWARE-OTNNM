function [plan, totalCost] = multi_obj_ot(marginal1, marginal2, costMatrix1, costMatrix2, gamma, lambda, max_iter)
% Multi-Objective Optimal Transport
% This function solves a multi-objective optimal transport problem using the Sinkhorn algorithm
% 
% Inputs:
%   marginal1, marginal2: marginal distributions (column vectors)
%   costMatrix1, costMatrix2: cost matrices for two objectives
%   gamma: regularization parameter
%   lambda: scalarization parameter (0 <= lambda <= 1)
%   max_iter: maximum number of iterations
%
% Outputs:
%   plan: optimal transport plan
%   totalCost: objective values for both cost matrices

% Initialize variables
n = length(marginal1);
m = length(marginal2);
u = ones(n, 1);  % Initialize u vector
v = ones(m, 1);  % Initialize v vector

% Compute Gibbs kernels
K1 = exp(-costMatrix1 / gamma);
K2 = exp(-costMatrix2 / gamma);
% Combine kernels using scalarization
K = lambda * K1 + (1-lambda) * K2;
% Alternative combination method (commented out):
% K = K1.^lambda .* K2.^(1-lambda);

% Sinkhorn iterations
for i = 1:max_iter
    u_prev = u;
    v_prev = v;
    
    % Update u
    u = marginal1 ./ (K * v);
    
    % Update v
    v = marginal2 ./ (K' * u);
    
    % Check convergence
    if norm(u./u_prev - 1, inf) < 1e-6 && norm(v./v_prev - 1, inf) < 1e-6
        break;
    end
end

% Compute optimal transport plan
plan = diag(u) * K * diag(v);

% Compute objective values for both cost matrices
totalCost = [sum(sum(plan .* costMatrix1)), sum(sum(plan .* costMatrix2))];
end