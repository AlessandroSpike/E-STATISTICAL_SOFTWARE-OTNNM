% Network Model Analysis and Visualization Script
%
% This script loads and analyzes results from different network models,
% allowing user selection of specific parameter combinations and
% visualizing the relationships between real and permuted network metrics.

clc;clearvars;close all

%% Load results based on user selection
list = {'weighted random graph','strength driven attachment',...
    'fitness based model directed'};
[indx,tf] = listdlg('ListString',list,'SelectionMode','single');

if indx==1
    load weighted_random_graph_results.mat
    prompt = {'Number of Nodes:';'Probability of Edge Formation:';...
        'Entropy:';'Number of Samples:'};
elseif indx==2
    load strength_driven_attachment_results.mat
     prompt = {'Number of Nodes:';'Number of Link for Each Node:';...
        'Entropy:';'Number of Samples:'};
elseif indx==3
    load fitness_based_model_directed_results.mat
     prompt = {'Number of Nodes:';'Probability of Edge Formation:';...
        'Entropy:';'Number of Samples:'};
else
    disp('You have not chosen.')
    disp('Weighted random graph selected by default')
    load weighted_random_graph_results.mat
    prompt = {'Number of Nodes:';'Probability of Edge Formation:';...
        'Entropy:';'Number of Samples:'};
end

%% Extract single param results for scatter plotting 
name = 'Select Parameters';

% Set up formats for parameter selection dialog
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

defaultanswer = {1,1,1,1};

% Display dialog for parameter selection
[answer, canceled] = inputsdlg(prompt, name, formats, defaultanswer);

% Extract selected parameter combination
param = [ParamCombinations(answer{1},1),...
    ParamCombinations(answer{2},2),ParamCombinations(answer{3},3)...
    ParamCombinations(answer{4},4)];
c = find(any(bsxfun(@minus,ParamCombinations,param),2)==0);
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

% Create figure for scatter plots
figure('Position', [100, 100, 800, 600])  % Increase figure size

% Subplot 1: In-Degree
subplot(2,2,1)
scatter(mean(indegE,2), indeg, 'filled', 'MarkerFaceAlpha', 0.5)
xlabel('Permuted')
ylabel('Real')
grid on
axis tight
title('In-Degree')
r = corr(mean(indegE,2), indeg);
annotation('textbox', 'String', sprintf('Correlation: %.3f', r), ...
           'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
           'Units', 'normalized', 'Position', [0.14, 0.82, 0.001, 0.001], 'FontSize', 10, 'Margin', 2)

% Subplot 2: Out-Degree
subplot(2,2,2)
scatter(mean(outdegE,2), outdeg, 'filled', 'MarkerFaceAlpha', 0.5)
xlabel('Permuted')
ylabel('Real')
grid on
axis tight
title('Out-Degree')
r = corr(mean(outdegE,2), outdeg);
annotation('textbox', 'String', sprintf('Correlation: %.3f', r), ...
           'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
           'Units', 'normalized', 'Position', [0.59, 0.82, 0.001, 0.001], 'FontSize', 10, 'Margin', 2)

% Subplot 3: In-Strength
subplot(2,2,3)
scatter(mean(instrE,2), instr, 'filled', 'MarkerFaceAlpha', 0.5)
xlabel('Permuted')
ylabel('Real')
grid on
axis tight
title('In-Strength')
r = corr(mean(instrE,2), instr);
annotation('textbox', 'String', sprintf('Correlation: %.3f', r), ...
           'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
           'Units', 'normalized', 'Position', [0.14, 0.37, 0.001, 0.001], 'FontSize', 10, 'Margin', 2)

% Subplot 4: Out-Strength
subplot(2,2,4)
scatter(mean(outstrE,2), outstr, 'filled', 'MarkerFaceAlpha', 0.5)
xlabel('Permuted')
ylabel('Real')
grid on
axis tight
title('Out-Strength')
r = corr(mean(outstrE,2), outstr);
annotation('textbox', 'String', sprintf('Correlation: %.3f', r), ...
           'FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
           'Units', 'normalized', 'Position', [0.59, 0.37, 0.001, 0.001], 'FontSize', 10, 'Margin', 2)

% Adjust subplot spacing
set(gcf, 'color', 'w');  % Set background color to white

%% Initialize containers for correlations
IndegCorr = zeros(size(ParamCombinations,1),1);
OutdegCorr = zeros(size(ParamCombinations,1),1);
InstrCorr = zeros(size(ParamCombinations,1),1);
OutstrCorr = zeros(size(ParamCombinations,1),1);

%% Compute correlations for all parameter combinations
for x = 1:size(ParamCombinations,1)
    IndegCorr(x) = corr(mean(PermIndeg{x},2), RealIndeg{x});
    OutdegCorr(x) = corr(mean(PermOutdeg{x},2), RealOutdeg{x});
    InstrCorr(x) = corr(mean(PermInstr{x},2), RealInstr{x});
    OutstrCorr(x) = corr(mean(PermOutstr{x},2), RealOutstr{x});
end

%% Find unique parameter combinations
combo = cell(size(ParamCombinations,2),1);
for x = 1:size(ParamCombinations,2)
    combo{x} = unique(ParamCombinations(:,x));
end

%% Plot effect of parameters on network metrics
figure('Position', [100, 100, 800, 600]);  % Adjust size as needed

for x = 1:size(ParamCombinations,2)
    subplot(2, 2, x)  % Assuming there are 4 parameters; adjust if different
    
    colors = {'#0072BD', '#D95319', '#EDB120', '#7E2F8E'};  % Define colors
    
    for y = 1:size(combo{x},1)
        aus = find(ParamCombinations(:,x)==combo{x}(y));
        IndegCorrParam = IndegCorr(aus,:);
        OutdegCorrParam = OutdegCorr(aus,:);
        InstrCorrParam = InstrCorr(aus,:);
        OutstrCorrParam = OutstrCorr(aus,:);
        
        % Reshape data for boxchart
        X = [IndegCorrParam; OutdegCorrParam; InstrCorrParam; OutstrCorrParam];
        G = repelem(1:4, length(IndegCorrParam));
        P = repmat(y, size(G));  % Group for different parameter values
        
        b = boxchart(G, X, 'GroupByColor', P);
        
        % Set colors for each group
        for i = 1:length(b)
            b(i).BoxFaceColor = colors{y};
            b(i).WhiskerLineColor = colors{y};
            b(i).MarkerColor = colors{y};
        end
        hold on
    end
    
    axis tight
    grid on
    title([char(ParamNames{x})])
    ylabel("Correlation")
    xticks(1:4)
    xticklabels({'In-Deg.','Out-Deg.','In-Str.','Out-Str.'})
    xlabel('Metrics')
    
    % Add legend
    legend(cellstr(num2str(combo{x})), 'Location', 'bestoutside')
end

sgtitle('Effect of Parameters on Network Metrics')