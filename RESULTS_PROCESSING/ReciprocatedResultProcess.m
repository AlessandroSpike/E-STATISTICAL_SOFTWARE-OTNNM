% Network Analysis and Visualization Script
%
% This script analyzes and visualizes the results of a network study focusing on 
% reciprocated strength-driven attachment. It performs the following tasks:
%
% 1. Loads pre-computed results from a .mat file
% 2. Allows user to select specific parameter combinations for detailed analysis
% 3. Creates scatter plots comparing real and permuted network metrics
% 4. Computes correlations between real and permuted metrics for all parameter combinations
% 5. Visualizes the effect of different parameters on network metrics using boxplots
%
% The script provides insights into how various parameters affect network properties
% such as in-degree, out-degree, in-strength, out-strength, and reciprocated links.

clc;clearvars;close all

%% Load results
% Load previously saved data from 'reciprocated_strength_driven_attachment_results.mat'
load reciprocated_strength_driven_attachment_results.mat

%% Extract single param results for scatter plotting 
% Set up dialog box for parameter selection
name = 'Select Parameters';
prompt = {'Number of Nodes:';'Number of Link for Each Node:';...
        'Entropy:';'Number of Samples:';'Cost weight:'};

% Define format for each parameter selection
formats(1,1).type   = 'list';
formats(1,1).style  = 'popupmenu';
formats(1,1).items  = unique(ParamCombinations(:,1));

formats(2,1).type   = 'list';
formats(2,1).style  = 'popupmenu';
formats(2,1).items  = unique(ParamCombinations(:,2));

formats(3,1).type   = 'list';
formats(3,1).style  = 'popupmenu';
formats(3,1).items  = unique(ParamCombinations(:,3));

formats(4,1).type   = 'list';
formats(4,1).style  = 'popupmenu';
formats(4,1).items  = unique(ParamCombinations(:,4));

formats(5,1).type   = 'list';
formats(5,1).style  = 'popupmenu';
formats(5,1).items  = unique(ParamCombinations(:,5));

defaultanswer = {1,1,1,1,1};

% Display dialog box and get user input
[answer, canceled] = inputsdlg(prompt, name, formats, defaultanswer);

% Extract selected parameter combination
param = [ParamCombinations(answer{1},1),...
    ParamCombinations(answer{2},2),ParamCombinations(answer{3},3)...
    ParamCombinations(answer{4},4),ParamCombinations(answer{5},5)];
c = find(any(bsxfun(@minus,ParamCombinations,param),2)==0);

% Display selected parameters
paramtable = [string(ParamNames)';string(param)];
disp(paramtable)

% Extract data for selected parameter combination
indegE = PermIndeg{c};
outdegE = PermOutdeg{c};
instrE = PermInstr{c};
outstrE = PermOutstr{c};
indeg = RealIndeg{c};
outdeg = RealOutdeg{c};
instr = RealInstr{c};
outstr = RealOutstr{c};
recdegE = PermReciproc{c};
degRec = RealReciproc{c};

% Create figure for scatter plots
figure('Position', [100, 100, 800, 600]) % Increase figure size

% Subplot 1: In-Degree
subplot(3,2,1)
scatter(mean(indegE,2), indeg, 'filled', 'MarkerFaceAlpha', 0.5)
xlabel('Permuted')
ylabel('Real')
grid on
axis tight
title('In-Degree')
r = corr(mean(indegE,2), indeg);
annotation('textbox', 'String', sprintf('Correlation: %.3f', r), ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
    'Units', 'normalized', 'Position', [0.14, 0.95, 0.001, 0.001], 'FontSize', 10, 'Margin', 2)

% Subplot 2: Out-Degree
subplot(3,2,2)
scatter(mean(outdegE,2), outdeg, 'filled', 'MarkerFaceAlpha', 0.5)
xlabel('Permuted')
ylabel('Real')
grid on
axis tight
title('Out-Degree')
r = corr(mean(outdegE,2), outdeg);
annotation('textbox', 'String', sprintf('Correlation: %.3f', r), ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
    'Units', 'normalized', 'Position', [0.59, 0.95, 0.001, 0.001], 'FontSize', 10, 'Margin', 2)

% Subplot 3: In-Strength
subplot(3,2,3)
scatter(mean(instrE,2), instr, 'filled', 'MarkerFaceAlpha', 0.5)
xlabel('Permuted')
ylabel('Real')
grid on
axis tight
title('In-Strength')
r = corr(mean(instrE,2), instr);
annotation('textbox', 'String', sprintf('Correlation: %.3f', r), ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
    'Units', 'normalized', 'Position', [0.14, 0.62, 0.001, 0.001], 'FontSize', 10, 'Margin', 2)

% Subplot 4: Out-Strength
subplot(3,2,4)
scatter(mean(outstrE,2), outstr, 'filled', 'MarkerFaceAlpha', 0.5)
xlabel('Permuted')
ylabel('Real')
grid on
axis tight
title('Out-Strength')
r = corr(mean(outstrE,2), outstr);
annotation('textbox', 'String', sprintf('Correlation: %.3f', r), ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
    'Units', 'normalized', 'Position', [0.59, 0.62, 0.001, 0.001], 'FontSize', 10, 'Margin', 2)

% Subplot 5: Reciprocated-Link
subplot(3,2,5)
scatter(mean(recdegE,2), degRec, 'filled', 'MarkerFaceAlpha', 0.5)
xlabel('Permuted')
ylabel('Real')
grid on
axis tight
title('Reciprocated-Link')
r = corr(mean(recdegE,2), degRec);
annotation('textbox', 'String', sprintf('Correlation: %.3f', r), ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
    'Units', 'normalized', 'Position', [0.14, 0.29, 0.001, 0.001], 'FontSize', 10, 'Margin', 2)

% Adjust subplot spacing
set(gcf, 'color', 'w'); % Set background color to white

%% Initialize containers for correlations
IndegCorr = zeros(size(ParamCombinations,1),1);
OutdegCorr = zeros(size(ParamCombinations,1),1);
InstrCorr = zeros(size(ParamCombinations,1),1);
OutstrCorr = zeros(size(ParamCombinations,1),1);
RecipCorr = zeros(size(ParamCombinations,1),1);

%% Compute correlations for all parameter combinations
for x = 1:size(ParamCombinations,1)
    IndegCorr(x) = corr(mean(PermIndeg{x},2),RealIndeg{x});
    OutdegCorr(x) = corr(mean(PermOutdeg{x},2),RealOutdeg{x});
    InstrCorr(x) = corr(mean(PermInstr{x},2),RealInstr{x});
    OutstrCorr(x) = corr(mean(PermOutstr{x},2),RealOutstr{x});
    RecipCorr(x) = corr(mean(PermReciproc{x},2),RealReciproc{x});
end

%% Find unique parameter combinations
combo = cell(size(ParamCombinations,2),1);
for x = 1:size(ParamCombinations,2)
    combo{x} = unique(ParamCombinations(:,x));
end

%% Plot effect of parameters on network metrics
figure('Position', [100, 100, 800, 600]);  % Adjust size as needed

for x = 1:size(ParamCombinations,2)
    subplot(2, 3, x)  % Create subplot for each parameter
    
    colors = {'#0072BD', '#D95319', '#EDB120', '#7E2F8E'};  % Define colors for different parameter values
    
    for y = 1:size(combo{x},1)
        % Find indices for current parameter value
        aus = find(ParamCombinations(:,x)==combo{x}(y));
        
        % Extract correlation data for current parameter value
        IndegCorrParam = IndegCorr(aus,:);
        OutdegCorrParam = OutdegCorr(aus,:);
        InstrCorrParam = InstrCorr(aus,:);
        OutstrCorrParam = OutstrCorr(aus,:);
        RecipCorrParam = RecipCorr(aus,:);
        
        % Reshape data for boxchart
        X = [IndegCorrParam; OutdegCorrParam; InstrCorrParam; OutstrCorrParam; RecipCorrParam];
        G = repelem(1:5, length(IndegCorrParam));
        P = repmat(y, size(G));  % Group for different parameter values
        
        % Create boxchart
        b = boxchart(G, X, 'GroupByColor', P);
        
        % Set colors for each group
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
    xticks(1:5)
    xticklabels({'In-Deg.','Out-Deg.','In-Str.','Out-Str.','Recip-Link.'})
    xlabel('Metrics')
    
    % Add legend for parameter values
    legend(cellstr(num2str(combo{x})), 'Location', 'bestoutside')
end

% Add overall title to the figure
sgtitle('Effect of Parameters on Network Metrics')