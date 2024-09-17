function [Ci, Q] = OTmodularity_dir(W, eta, plan)
    % Optimal Transport Modularity for Directed Networks
    % This function finds the optimal community structure and calculates modularity for a directed network
    % Inputs:
    %   W: directed weighted/binary connection matrix
    %   eta: resolution parameter (optional, default = 1)
    %   plan: optimal transport plan
    % Outputs:
    %   Ci: optimal community structure
    %   Q: maximized modularity

    % Set default eta if not provided
    if ~exist('eta', 'var')
        eta = 1;
    end

    % Initialize variables
    N = length(W);
    Ki = sum(W, 1);  % in-degree
    m = sum(Ki);

    % Compute modularity matrix
    b = W - eta * plan;
    B = b + b.';  % Symmetrize the modularity matrix

    % Initialize community structure
    Ci = ones(N, 1);
    cn = 1;  % Number of communities
    U = [1 0];  % Array of unexamined communities

    ind = 1:N;
    Bg = B;
    Ng = N;

    % Main community detection loop
    while U(1)
        % Compute eigenvector corresponding to largest positive eigenvalue
        [V, D] = eig(Bg);
        [~, i1] = max(real(diag(D)));
        v1 = V(:, i1);

        % Initial community division based on eigenvector sign
        S = ones(Ng, 1);
        S(v1 < 0) = -1;
        q = S.' * Bg * S;  % Contribution to modularity

        if q > 1e-10  % Community is divisible
            qmax = q;
            Bg(logical(eye(Ng))) = 0;  % Modify Bg for fine-tuning
            indg = ones(Ng, 1);
            Sit = S;

            % Fine-tuning loop
            while any(indg)
                Qit = qmax - 4 * Sit .* (Bg * Sit);
                [qmax, imax] = max(Qit .* indg);
                Sit(imax) = -Sit(imax);
                indg(imax) = nan;
                if qmax > q
                    q = qmax;
                    S = Sit;
                end
            end

            % Update community structure
            if abs(sum(S)) == Ng  % Unsuccessful splitting
                U(1) = [];
            else
                cn = cn + 1;
                Ci(ind(S == 1)) = U(1);
                Ci(ind(S == -1)) = cn;
                U = [cn U];
            end
        else  % Community is indivisible
            U(1) = [];
        end

        % Prepare for next iteration
        ind = find(Ci == U(1));
        bg = B(ind, ind);
        Bg = bg - diag(sum(bg));
        Ng = length(ind);
    end

    % Compute final modularity
    s = Ci(:, ones(1, N));
    Q = ~(s - s.') .* B / (2 * m);
    Q = sum(Q(:));
end