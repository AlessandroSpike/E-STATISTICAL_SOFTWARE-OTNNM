function visualize_network(W, options, community_labels)
% VISUALIZE_NETWORK Visualize a network given its adjacency matrix
%   W: Adjacency matrix (can be binary or weighted)
%   options: (Optional) Structure with visualization options
%   community_labels: (Optional) Vector of community labels for nodes
%
% Example usage:
% 1. With only adjacency matrix:
%    visualize_network(W)
% 
% 2. With custom options:
%    options = struct('Layout', 'force', 'EdgeWidthScale', 2, 'NodeSizeScale', 20, 'ColorMap', 'jet', 'Directed', true);
%    visualize_network(W, options)
%
% 3. With community labels:
%    community_labels = [1, 1, 2, 2, 3, 3];
%    visualize_network(W, options, community_labels)

    % Check inputs and set defaults
    if nargin < 2
        options = struct();
    end
    if nargin < 3
        community_labels = [];
    end

    % Set default options
    options = set_default_options(options);

    % Create graph object
    if options.Directed
        G = digraph(W);
    else
        G = graph(W);
    end

    % Compute node sizes
    node_sizes = compute_node_sizes(W, options);

    % Compute edge widths and colors
    [edge_widths, edge_colors] = compute_edge_properties(G, options);

    % Determine node colors
    if ~isempty(community_labels)
        node_colors = 10*community_labels;
        color_data = 'community';
    else
        node_colors = sum(abs(W), 2);  % Node strength
        color_data = 'strength';
    end

    % Plot the graph
    figure;
    h = plot(G, 'Layout', options.Layout);

    % Set node properties
    h.NodeCData = node_colors;
    h.MarkerSize = node_sizes;

    % Set edge properties
    h.LineWidth = edge_widths;
    h.EdgeCData = edge_colors;

    % Set axis properties
    axis off;
    title('Network Visualization');

    % Add colorbar for edges
    colormap(options.ColorMap);
    c_edge = colorbar('Location', 'southoutside');
    c_edge.Label.String = 'Edge Weight';

    % Add colorbar for nodes if using strength-based coloring
    if strcmp(color_data, 'strength')
        c_node = colorbar('Location', 'eastoutside');
        c_node.Label.String = 'Node Strength';
    end

    
end

function options = set_default_options(options)
    % Set default options if not provided
    default_options = struct('Layout', 'force', ...
                             'EdgeWidthScale', 2, ...
                             'NodeSizeScale', 20, ...
                             'ColorMap', 'hot', ...
                             'Directed', false);
    options = setstructfields(options,default_options);
end

function node_sizes = compute_node_sizes(W, options)
    node_strengths = sum(abs(W), 2);  % Use absolute values for negative weights
    node_sizes = normalize_sizes(node_strengths, options.NodeSizeScale);
end

function [edge_widths, edge_colors] = compute_edge_properties(G, options)
    edge_weights = abs(G.Edges.Weight);  % Use absolute values for negative weights
    edge_widths = normalize_sizes(edge_weights, options.EdgeWidthScale);
    edge_colors = edge_weights;
end

function normalized = normalize_sizes(sizes, scale)
    % Normalize sizes to a range suitable for visualization
    min_size = 1;
    normalized = (sizes - min(sizes)) / (max(sizes) - min(sizes));
    normalized = normalized * scale + min_size;
end

function S = setstructfields(S, T)
    % Set fields in S from T, without overwriting existing fields
    names = fieldnames(T);
    for i = 1:length(names)
        if ~isfield(S, names{i})
            S.(names{i}) = T.(names{i});
        end
    end
end