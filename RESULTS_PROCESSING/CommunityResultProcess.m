% Network Parameter Analysis Script
% 
% This script analyzes and visualizes the effects of various parameters on 
% network metrics in a complex network study. It performs the following tasks:
%
% 1. Loads pre-computed results from a .mat file
% 2. Identifies unique parameter combinations used in the study
% 3. Creates a series of boxplots to visualize how different parameters 
%    affect network metrics such as modularity and variation
% 4. Each subplot represents a different parameter, showing how its values 
%    influence the distribution of metrics
% 5. The metrics visualized include:
%    - Modularity for OT null model (OT) and Newmann traditional modularity
%    - Variation of information for OT null model (OT) and Newmann traditional modularity
%
% The resulting figure provides a comprehensive view of how each parameter 
% impacts network structure and behavior across different metrics and models.

clc;clearvars;close all

%% Load results
% Load previously saved data from 'community_results.mat'
load community_results.mat

%% Find unique parameter combinations
combo = cell(size(ParamCombinations,2),1);
for x = 1:size(ParamCombinations,2)
    % Extract unique values for each parameter
    combo{x} = unique(ParamCombinations(:,x));
end

%% Plot results
% Create a new figure with specified position and size
figure('Position', [100, 100, 800, 600]);  % Adjust size as needed

% Loop through each parameter
for x = 1:size(ParamCombinations,2)
    % Create subplot for each parameter
    subplot(2, 3, x) 
    
    % Define colors for different parameter values
    colors = {'#0072BD', '#D95319', '#EDB120', '#7E2F8E'};  
    
    % Loop through each unique value of the current parameter
    for y = 1:size(combo{x},1)
        % Find indices where the parameter has this value
        aus = find(ParamCombinations(:,x) == combo{x}(y));

        % Extract and reshape data for plotting
        ModulOT = PermModul(aus,:);
        ModulOT = ModulOT(:);
        ModulNeu = RealModul(aus,:);
        ModulNeu = ModulNeu(:);
        VariationOT = PermVariation(aus,:);
        VariationOT = VariationOT(:);
        VariationNeu = RealVariation(aus,:);
        VariationNeu = VariationNeu(:);
        
        % Prepare data for boxchart
        X = [ModulOT;ModulNeu;VariationOT;VariationNeu];
        G = repelem(1:4, length(ModulOT));
        P = repmat(y, size(G));  % Group for different parameter values
        
        % Create boxchart
        b = boxchart(G, X, 'GroupByColor', P);
        
        % Set colors for each group in the boxchart
        for i = 1:length(b)
            b(i).BoxFaceColor = colors{y};
            b(i).WhiskerLineColor = colors{y};
            b(i).MarkerColor = colors{y};
        end
        hold on
    end
    
    % Adjust plot appearance
    axis tight
    grid on
    title([char(ParamNames{x})])
    ylabel("Correlation")
    xticks(1:6)
    xticklabels({'Modul.-OT','Modul.-New.','Variation-OT','Variation-New.'})
    xlabel('Metrics')
    
    % Add legend for parameter values
    legend(cellstr(num2str(combo{x})), 'Location', 'bestoutside')
   
end

% Add overall title to the figure
sgtitle('Effect of Parameters on Network Metrics')