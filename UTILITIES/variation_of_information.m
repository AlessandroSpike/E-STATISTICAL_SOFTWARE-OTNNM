function [VIn] = variation_of_information(Cx, Cy)
% VARIATION_OF_INFORMATION Calculates the normalized variation of information between community partitions
%
% This function quantifies information-theoretic distance (normalized
% variation of information) between community partitions.
%
% Inputs:
%   Cx - Community partition vector or matrix of n rows and p columns,
%        where n is the number of network nodes, and p is the number of input
%        community partitions (in the case of vector input p=1).
%   Cy - (Optional) Community partition vector or matrix of n rows and q columns.
%        If omitted, the partition distance is computed between all pairwise partitions of Cx.
%
% Outputs:
%   VIn - Normalized variation of information ([p, q] matrix)
%
% Notes:
%   VIn = [H(X) + H(Y) - 2MI(X, Y)]/log(n)
%   where H is the entropy and MI is the mutual information

% Check if only one input is provided
s = (nargin==1);
if s
    Cy = Cx;
    d = 10.^ceil(log10(double(1 + max( Cx(:)) )));
else
    d = 10.^ceil(log10(double(1 + max([Cx(:);Cy(:)]) )));
end

% Check if inputs are positive integers
if ~isequal([Cx(:);Cy(:)], int64([Cx(:);Cy(:)])) || min([Cx(:);Cy(:)])<=0
    error('Input partitions must contain only positive integers.')
end

% Get dimensions of input
[n, p] = size(Cx);

% Calculate entropy for Cx
HX = zeros(p, 1);
for i = 1:p
    Px = nonzeros(accumarray(Cx(:, i), 1)) / n;                     % P(x)
    HX(i) = - sum(Px .* log(Px));                                   % H(x)
end

% Calculate entropy for Cy (or use HX if Cy is not provided)
if s
    q = p;
    HY = HX;
else
    [n_, q] = size(Cy);
    assert(n == n_);
    HY = zeros(q, 1);
    for j = 1:q
        Py = nonzeros(accumarray(Cy(:, j), 1)) / n;                 % P(y)
        HY(j) = - sum(Py .* log(Py));                               % H(y)
    end
end

% Initialize output matrix
VIn = zeros(p, q);

% Calculate normalized variation of information
for i = 1:p
    j_idx = (s * (i - 1) + 1):q;
    for j = j_idx
        Pxy = nonzeros(accumarray(d*Cx(:, i) + Cy(:, j), 1)) / n;       % P(x,y)
        Hxy = -sum(Pxy .* log(Pxy));                                % H(x,y)
        VIn(i, j) = (2 * Hxy - HX(i) - HY(j)) / log(n);             % VIn
    end
    % If only one input, fill in symmetric part of matrix
    if s
        VIn(j_idx, i) = VIn(i, j_idx);
    end
end