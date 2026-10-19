function plotascascores(factormodel, sccomp, ascamodel)
% plotascascores Plots ASCA scores plot
%  The function plots the values of ASCA scores eitehr as a bar plot (if 
%  model includes only one component, or as 2D scatterplot if at least two 
%  components are extracted. 
%  First plot plots only the scores for each factor level, while the second 
%  one plots the scores after projection of the residuals).
%  
%
%  INPUTS:
%   factormodel = field of ASCA model corresponding to the design term
%  OPTIONAL INPUTS:
%        sccomp    = Vector with the number of components to be plotted (if
%                    empty 1:min(ncomp,2)).
%        ascamodel = full ASCA model. If provided, the function will build a
%                    legend with factor-level labels from the ORIGINAL
%                    design input stored in the model (e.g., model.DesignInput
%                    and model.DesignInfo).
%  I/O:
%       plotascascores(factormodel, sccomp);
%       plotascascores(factormodel);
%       plotascascores(factormodel, [], model);
%
%
% Author: Federico Marini
% Version: 15/05/2020

if nargin < 2 || isempty(sccomp) || max(sccomp) > size(factormodel.SCA.Model.Scores,2)
    sccomp=1:min(size(factormodel.SCA.Model.Scores,2), 2);
end

% If the full model is provided, build legend labels based on the original
% factor levels for the factors involved in this term.
haveLegend = (nargin >= 3) && ~isempty(ascamodel) && isstruct(ascamodel);
legendLabels = {};
if haveLegend
    try
        legendLabels = local_build_legend_labels(factormodel, ascamodel);
        haveLegend = ~isempty(legendLabels);
    catch
        haveLegend = false;
    end
end


[~,~,lindex]=unique(factormodel.DesignMatrix,'rows', 'stable');
Mark={'o','s','d','^','v','>','<','p','h','+','*','.','x'};
colr=lines(7);

if length(sccomp)>1
    figure
    hold on

    h = gobjects(max(lindex),1);
    
    for i=1:max(lindex)
        f1=rem(i,13); 
        if f1==0 
            f1=13;
        end
        f2=rem(i,7); 
        if f2==0 
            f2=7;
        end
        
        h(i)=plot(factormodel.SCA.Model.Scores(lindex==i,sccomp(1)), factormodel.SCA.Model.Scores(lindex==i,sccomp(2)),Mark{f1}, 'markerfacecolor', colr(f2,:), 'markeredgecolor', colr(f2,:), 'markersize', 10);
    end
    
    xx=axis; 
    plot(xx(1:2), [0 0], '--k', 'linewidth', 1.2)
    plot([0 0],xx(3:4),  '--k', 'linewidth', 1.2)
    xlabel(['Scores on SC', num2str(sccomp(1))], 'fontsize', 12, 'fontweight', 'bold')
    ylabel(['Scores on SC', num2str(sccomp(2))], 'fontsize', 12, 'fontweight', 'bold')
    set(gca, 'linewidth', 1.5, 'fontsize', 12, 'fontweight', 'bold')

    if haveLegend
        legend(h, legendLabels, 'Interpreter','none', 'Location','best');
    end
    
    
    figure
    hold on

    h = gobjects(max(lindex),1);
    
    for i=1:max(lindex)
        f1=rem(i,13); 
        if f1==0 
            f1=13;
        end
        f2=rem(i,7); 
        if f2==0 
            f2=7;
        end
        
        h(i)=plot(factormodel.SCA.Model.Scores(lindex==i,sccomp(1)), factormodel.SCA.Model.Scores(lindex==i,sccomp(2)),Mark{f1}, 'markerfacecolor', colr(f2,:), 'markeredgecolor', colr(f2,:), 'markersize', 10);
        plot(factormodel.SCA.Model.ScoreswithRes(lindex==i,sccomp(1)), factormodel.SCA.Model.ScoreswithRes(lindex==i,sccomp(2)),Mark{f1}, 'markeredgecolor', colr(f2,:), 'markersize', 6)
    end
    
    xx=axis; 
    plot(xx(1:2), [0 0], '--k', 'linewidth', 1.2)
    plot([0 0],xx(3:4),  '--k', 'linewidth', 1.2)
    xlabel(['Scores on SC', num2str(sccomp(1))], 'fontsize', 12, 'fontweight', 'bold')
    ylabel(['Scores on SC', num2str(sccomp(2))], 'fontsize', 12, 'fontweight', 'bold')
    set(gca, 'linewidth', 1.5, 'fontsize', 12, 'fontweight', 'bold')

    if haveLegend
        legend(h, legendLabels, 'Interpreter','none', 'Location','best');
    end
    
    
else
    figure
    hold on

    h = gobjects(max(lindex),1);
    for i=1:max(lindex)
        f2=rem(i,7); 
        if f2==0 
            f2=7;
        end
        
        h(i)=bar(find(lindex==i), factormodel.SCA.Model.Scores(lindex==i,sccomp(1)),'facecolor', colr(f2,:), 'edgecolor', colr(f2,:));
        
    
    end
    
    axis tight
    xlabel('Sample Index', 'fontsize', 12, 'fontweight', 'bold')
    ylabel(['Scores on SC', num2str(sccomp(1))], 'fontsize', 12, 'fontweight', 'bold')
    set(gca, 'linewidth', 1.5, 'fontsize', 12, 'fontweight', 'bold')

    if haveLegend
        legend(h, legendLabels, 'Interpreter','none', 'Location','best');
    end
    
    figure
    hold on

    h = gobjects(max(lindex),1);
    for i=1:max(lindex)
        f2=rem(i,7); 
        if f2==0 
            f2=7;
        end
        
        h(i)=bar(find(lindex==i), factormodel.SCA.Model.ScoreswithRes(lindex==i,sccomp(1)),'facecolor', colr(f2,:), 'edgecolor', colr(f2,:));

        
    end
    
    axis tight
    xlabel('Sample Index', 'fontsize', 12, 'fontweight', 'bold')
    ylabel(['Scores on SC', num2str(sccomp(1))], 'fontsize', 12, 'fontweight', 'bold')
    set(gca, 'linewidth', 1.5, 'fontsize', 12, 'fontweight', 'bold')

    if haveLegend
        legend(h, legendLabels, 'Interpreter','none', 'Location','best');
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function labels = local_build_legend_labels(factormodel, model)
% Build legend labels from the ORIGINAL design input stored in the model.
% The grouping is inferred from unique rows of factormodel.DesignMatrix.
% For each group, we take the first sample and report the factor levels for
% the factors involved in this term.

labels = {};

% Need the original design input (numeric matrix, cell array, or string array)
if isfield(model,'DesignInput')
    D0 = model.DesignInput;
elseif isfield(model,'Design')
    % fallback: coded design
    D0 = model.Design;
else
    return;
end

% Determine which factors are involved in this term.
% Prefer explicit term members if available.
terms = {};
if isfield(factormodel,'Members')
    terms = factormodel.Members;
elseif isfield(factormodel,'DesignTerms')
    % DesignTerms might be numeric indices of factors
    dt = factormodel.DesignTerms;
    if isnumeric(dt)
        terms = arrayfun(@(x) char(64+x), dt(:)', 'UniformOutput', false);
    end
end
if isempty(terms)
    % last resort: use TermLabel if present
    if isfield(factormodel,'EffectLabel')
        terms = {factormodel.EffectLabel};
    else
        return;
    end
end

% Union of base factor letters used by any member term
letters = '';
for k = 1:numel(terms)
    t = char(terms{k});
    t = t(isstrprop(t,'alpha'));
    letters = [letters t]; %#ok<AGROW>
end
letters = unique(letters, 'stable');
if isempty(letters)
    return;
end
fidx = double(letters) - 64;  % 'A'->1

% Grouping index used in the plot
[~,~,lindex] = unique(factormodel.DesignMatrix,'rows','stable');
ng = max(lindex);
labels = cell(ng,1);

for g = 1:ng
    ii = find(lindex==g, 1, 'first');
    if isempty(ii)
        labels{g} = sprintf('group %d', g);
        continue;
    end
    labels{g} = local_level_string(D0, ii, fidx, model);
end


function s = local_level_string(D0, rowIdx, fidx, model)
% Create a compact label for a given sample row based on factor indices.

parts = cell(1, numel(fidx));
for k = 1:numel(fidx)
    j = fidx(k);
    fname = char(64+j);

    % factor name (if present)
    if isfield(model,'DesignInfo') && isfield(model.DesignInfo,'Factors') && numel(model.DesignInfo.Factors) >= j
        if isfield(model.DesignInfo.Factors(j), 'Name') && ~isempty(model.DesignInfo.Factors(j).Name)
            fname = model.DesignInfo.Factors(j).Name;
        end
    end

    % extract original level value
    if isnumeric(D0)
        v = D0(rowIdx, j);
        lv = num2str(v);
    elseif iscell(D0)
        lv = D0{rowIdx, j};
        if isnumeric(lv)
            lv = num2str(lv);
        else
            lv = char(string(lv));
        end
    else
        % string array / categorical
        lv = char(string(D0(rowIdx, j)));
    end

    parts{k} = sprintf('%s=%s', fname, lv);
end
s = strjoin(parts, ', ');


        
        