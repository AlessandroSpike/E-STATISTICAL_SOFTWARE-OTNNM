function [in_strength, out_strength, in_degree, out_degree, rec_degree] = compute_network_properties(W)
% compute_network_properties Calculates network properties from adjacency matrix
%
% Inputs:
%   W: NxN adjacency matrix
%
% Outputs:
%   in_strength: Nx1 vector of in-strengths
%   out_strength: Nx1 vector of out-strengths
%   in_degree: Nx1 vector of in-degrees
%   out_degree: Nx1 vector of out-degrees
%   rec_degree: Nx1 vector of reciprocated links

% Ensure W is a 2D matrix
if ~ismatrix(W)
    error('Input A must be a 2D matrix');
end

% Compute in-strength and out-strength
in_strength = sum(W, 1)';  % sum along columns
out_strength = sum(W, 2);  % sum along rows

% Compute in-degree and out-degree
in_degree = sum(W > 0, 1)';  % sum of non-zero elements along columns
out_degree = sum(W > 0, 2);  % sum of non-zero elements along rows

% Compute reciprocated links for each node
R = W & W';
rec_degree = sum(double(R), 1)'; 

% Ensure all outputs are column vectors
in_strength = reshape(in_strength, [], 1);
out_strength = reshape(out_strength, [], 1);
in_degree = reshape(in_degree, [], 1);
out_degree = reshape(out_degree, [], 1);
rec_degree = reshape(rec_degree, [], 1);

end