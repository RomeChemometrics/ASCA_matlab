function plotascascores_checked(factormodel, sccomp, ascamodel)
% plotascascores_checked  Plot ASCA/SCA scores with correct cell grouping.
%
% This version fixes a common problem in ASCA score plots: grouping samples
% according to unique rows of factormodel.DesignMatrix is not always
% equivalent to grouping samples according to the real experimental cells.
% This is especially problematic for interactions when effect coding is used.
%
% INPUTS
%   factormodel : one term of an ASCA model, e.g. model.XA, model.XAB, ...
%   sccomp      : components to plot. If empty, 1:min(2,ncomp) is used.
%   ascamodel   : optional full ASCA model. Strongly recommended, because
%                 it contains the original design matrix. If provided, the
%                 grouping is based on the original factor levels defining
%                 the term, not on the contrast-coded design matrix.
%
% I/O
%   plotascascores_checked(model.XA, [], model)
%   plotascascores_checked(model.XAB, [1 2], model)
%
% Main correction with respect to the original function:
%   The plotting groups are built from the original design cells:
%       model.DesignInput(:, factormodel.DesignTerms)
%   rather than from:
%       factormodel.DesignMatrix
%   because effect-coded interaction matrices may assign the same coded row
%   to different cells or otherwise not correspond one-to-one to cells.
%
% Notes:
%   Figure 1 shows only the model term scores, averaged by cell.
%   Figure 2 shows the cell model score plus residual-projected sample scores.
%
% R2015b-oriented: avoids contains/isstring/string/gobjects.

% -------------------------------------------------------------------------
% Check SCA fields
% -------------------------------------------------------------------------
if ~isfield(factormodel,'SCA') || ~isfield(factormodel.SCA,'Model') || ...
        ~isfield(factormodel.SCA.Model,'Scores')
    error('plotascascores_checked:MissingScores', ...
        'factormodel must contain factormodel.SCA.Model.Scores.');
end

Scores = factormodel.SCA.Model.Scores;
ncomp  = size(Scores,2);

if nargin < 2 || isempty(sccomp) || max(sccomp) > ncomp
    sccomp = 1:min(ncomp,2);
end

if isempty(sccomp)
    error('plotascascores_checked:NoComponents', ...
        'No SCA components are available for this term.');
end

% Use only one or two components.
if length(sccomp) > 2
    sccomp = sccomp(1:2);
end

% -------------------------------------------------------------------------
% Build grouping based on original experimental cells whenever possible.
% -------------------------------------------------------------------------
if nargin >= 3 && ~isempty(ascamodel) && isstruct(ascamodel)
    [gindex, labels, groupName] = local_group_from_original_design(factormodel, ascamodel);
else
    [gindex, labels, groupName] = local_group_from_design_matrix(factormodel);
end

ng = max(gindex);

% -------------------------------------------------------------------------
% Compute one model score point per cell.
% For a correctly fitted ASCA model, Scores should normally be identical
% within the same cell for the effect model. We nevertheless average them
% to avoid plotting multiple superimposed "average" points.
% -------------------------------------------------------------------------
ScoreMean = zeros(ng, length(sccomp));
for g = 1:ng
    ScoreMean(g,:) = mean(Scores(gindex == g, sccomp), 1);
end

% Residual-projected scores, if available.
hasRes = isfield(factormodel.SCA.Model,'ScoreswithRes') && ...
         ~isempty(factormodel.SCA.Model.ScoreswithRes);
if hasRes
    ScoresRes = factormodel.SCA.Model.ScoreswithRes;
else
    ScoresRes = [];
end

% -------------------------------------------------------------------------
% Plot settings
% -------------------------------------------------------------------------
Mark = {'o','s','d','^','v','>','<','p','h','+','*','.','x'};
colr = lines(7);

% MATLAB R2015b compatible graphics handle initialization
h = zeros(ng,1);

% =========================================================================
% Case 1: two-dimensional score plot
% =========================================================================
if length(sccomp) > 1

    % ---------------------------------------------------------------------
    % Figure 1: one model score point per experimental cell
    % ---------------------------------------------------------------------
    figure('units', 'normalized', 'Position',[0.15 0.15 0.8 0.8])
    hold on

    for g = 1:ng
        [mk, cc] = local_marker_color(g, Mark, colr);
        h(g) = plot(ScoreMean(g,1), ScoreMean(g,2), mk, ...
            'markerfacecolor', cc, ...
            'markeredgecolor', cc, ...
            'markersize', 15, ...
            'linestyle', 'none');
    end

    local_add_axes
    xlabel(['Scores on SC', num2str(sccomp(1))], ...
        'fontsize', 16, 'fontweight', 'bold')
    ylabel(['Scores on SC', num2str(sccomp(2))], ...
        'fontsize', 16, 'fontweight', 'bold')
    title(['ASCA model scores grouped by ', groupName], ...
        'fontsize', 16, 'fontweight', 'bold')
    set(gca, 'linewidth', 1.5, 'fontsize', 16, 'fontweight', 'bold')

    if ~isempty(labels)
        legend(h, labels, 'Interpreter','none', 'Location','best');
    end

    % ---------------------------------------------------------------------
    % Figure 2: model cell score + residual-projected individual samples
    % ---------------------------------------------------------------------
    figure('units', 'normalized', 'Position',[0.15 0.15 0.8 0.8])
    
    hold on

    h = zeros(ng,1);

    for g = 1:ng
        [mk, cc] = local_marker_color(g, Mark, colr);

        % Residual-projected individual samples for this cell
        if hasRes
            plot(ScoresRes(gindex == g, sccomp(1)), ...
                 ScoresRes(gindex == g, sccomp(2)), mk, ...
                 'markeredgecolor', cc, ...
                 'markersize', 10, ...
                 'linestyle', 'none');
        end

        % Cell model score, plotted once
        h(g) = plot(ScoreMean(g,1), ScoreMean(g,2), mk, ...
            'markerfacecolor', cc, ...
            'markeredgecolor', cc, ...
            'markersize', 14, ...
            'linestyle', 'none');
    end

    local_add_axes
    xlabel(['Scores on SC', num2str(sccomp(1))], ...
        'fontsize', 16, 'fontweight', 'bold')
    ylabel(['Scores on SC', num2str(sccomp(2))], ...
        'fontsize', 16, 'fontweight', 'bold')
    title(['ASCA scores with residual projection grouped by ', groupName], ...
        'fontsize', 16, 'fontweight', 'bold')
    set(gca, 'linewidth', 1.5, 'fontsize', 16, 'fontweight', 'bold')

    if ~isempty(labels)
        legend(h, labels, 'Interpreter','none', 'Location','best');
    end

% =========================================================================
% Case 2: one-dimensional score plot
% =========================================================================
else

    % ---------------------------------------------------------------------
    % Figure 1: one bar per experimental cell
    % ---------------------------------------------------------------------
    figure
    hold on

    for g = 1:ng
        [~, cc] = local_marker_color(g, Mark, colr);
        h(g) = bar(g, ScoreMean(g,1), ...
            'facecolor', cc, ...
            'edgecolor', cc);
    end

    axis tight
    xlabel(['Cell index: ', groupName], ...
        'fontsize', 16, 'fontweight', 'bold')
    ylabel(['Scores on SC', num2str(sccomp(1))], ...
        'fontsize', 16, 'fontweight', 'bold')
    title('ASCA model scores: one average point/bar per cell', ...
        'fontsize', 16, 'fontweight', 'bold')
    set(gca, 'linewidth', 1.5, 'fontsize', 16, 'fontweight', 'bold')

    if ~isempty(labels)
        legend(h, labels, 'Interpreter','none', 'Location','best');
    end

    % ---------------------------------------------------------------------
    % Figure 2: individual residual-projected scores plus cell model score
    % ---------------------------------------------------------------------
    figure
    hold on

    h = zeros(ng,1);

    for g = 1:ng
        [mk, cc] = local_marker_color(g, Mark, colr);

        % Plot residual-projected individual samples horizontally around
        % their cell index.
        if hasRes
            x = g * ones(sum(gindex == g),1);
            plot(x, ScoresRes(gindex == g, sccomp(1)), mk, ...
                'markeredgecolor', cc, ...
                'markersize', 6, ...
                'linestyle', 'none');
        end

        % Plot model score once as a filled marker
        h(g) = plot(g, ScoreMean(g,1), mk, ...
            'markerfacecolor', cc, ...
            'markeredgecolor', cc, ...
            'markersize', 11, ...
            'linestyle', 'none');
    end

    axis tight
    xlabel(['Cell index: ', groupName], ...
        'fontsize', 16, 'fontweight', 'bold')
    ylabel(['Scores on SC', num2str(sccomp(1))], ...
        'fontsize', 16, 'fontweight', 'bold')
    title('ASCA scores with residual projection', ...
        'fontsize', 16, 'fontweight', 'bold')
    set(gca, 'linewidth', 1.5, 'fontsize', 16, 'fontweight', 'bold')

    if ~isempty(labels)
        legend(h, labels, 'Interpreter','none', 'Location','best');
    end
end


% =========================================================================
% Local helper functions
% =========================================================================

function [gindex, labels, groupName] = local_group_from_original_design(factormodel, model)
% Build group index from the original design matrix and the factors involved
% in the selected ASCA term.

% Determine which original factors define the term.
fidx = [];

if isfield(factormodel,'DesignTerms') && isnumeric(factormodel.DesignTerms)
    fidx = factormodel.DesignTerms(:)';
end

% Fallback: derive factor indices from a label such as 'A', 'AB', 'A_B'.
if isempty(fidx) && isfield(factormodel,'EffectLabel')
    lab = factormodel.EffectLabel;
    fidx = local_factor_indices_from_label(lab);
end

if isempty(fidx)
    % If the factors cannot be determined, use the term design matrix.
    [gindex, labels, groupName] = local_group_from_design_matrix(factormodel);
    groupName = [groupName, ' (fallback)'];
    return
end

% Get the original design if possible; otherwise use the coded design.
if isfield(model,'DesignInput')
    D0 = model.DesignInput;
elseif isfield(model,'Design')
    D0 = model.Design;
else
    [gindex, labels, groupName] = local_group_from_design_matrix(factormodel);
    groupName = [groupName, ' (fallback)'];
    return
end

% Protect against invalid factor indices.
nf = size(D0,2);
fidx = fidx(fidx >= 1 & fidx <= nf);
if isempty(fidx)
    [gindex, labels, groupName] = local_group_from_design_matrix(factormodel);
    groupName = [groupName, ' (fallback)'];
    return
end

% Build grouping from the original design columns involved in the term.
if isnumeric(D0)
    Dcell = D0(:, fidx);
    [~,~,gindex] = unique(Dcell, 'rows', 'stable');
else
    Dcell = D0(:, fidx);
    [~,~,gindex] = local_unique_cell_rows(Dcell);
end

% Build readable labels.
ng = max(gindex);
labels = cell(ng,1);
for g = 1:ng
    ii = find(gindex == g, 1, 'first');
    labels{g} = local_level_string(D0, ii, fidx, model);
end

groupName = local_factor_name_string(fidx);


function [gindex, labels, groupName] = local_group_from_design_matrix(factormodel)
% Fallback grouping based on the ASCA term design matrix.
% This is not recommended for interactions with effect coding.

if isfield(factormodel,'DesignMatrix')
    [~,~,gindex] = unique(factormodel.DesignMatrix, 'rows', 'stable');
else
    ns = size(factormodel.SCA.Model.Scores,1);
    gindex = (1:ns)';
end

ng = max(gindex);
labels = cell(ng,1);
for g = 1:ng
    labels{g} = ['group ', num2str(g)];
end
groupName = 'term design matrix';


function fidx = local_factor_indices_from_label(lab)
% Convert labels such as 'A', 'AB', 'A_B', 'A+AB' into factor indices.

if isempty(lab)
    fidx = [];
    return
end

lab = char(lab);
letters = lab(isletter(lab));
letters = unique(letters, 'stable');
fidx = double(upper(letters)) - 64;


function s = local_level_string(D0, rowIdx, fidx, model)
% Create a label such as A=low, B=high for one experimental cell.

parts = cell(1, numel(fidx));

for k = 1:numel(fidx)
    j = fidx(k);
    fname = char(64+j);

    % Use stored factor name if available.
    if isfield(model,'DesignInfo') && isfield(model.DesignInfo,'Factors') && ...
            numel(model.DesignInfo.Factors) >= j && ...
            isfield(model.DesignInfo.Factors(j),'Name') && ...
            ~isempty(model.DesignInfo.Factors(j).Name)
        fname = model.DesignInfo.Factors(j).Name;
    end

    if isnumeric(D0)
        lv = num2str(D0(rowIdx,j));
    elseif iscell(D0)
        val = D0{rowIdx,j};
        if isnumeric(val)
            lv = num2str(val);
        else
            lv = char(val);
        end
    else
        % Basic fallback for other array-like inputs.
        val = D0(rowIdx,j);
        if isnumeric(val)
            lv = num2str(val);
        else
            lv = char(val);
        end
    end

    parts{k} = [fname, '=', lv];
end

s = local_strjoin(parts, ', ');


function s = local_factor_name_string(fidx)
% Return a compact name for the involved factors, e.g. A+B or A+B+C.

parts = cell(1, numel(fidx));
for k = 1:numel(fidx)
    parts{k} = char(64 + fidx(k));
end
s = local_strjoin(parts, '+');


function [C, ia, ic] = local_unique_cell_rows(A)
% R2015b-safe unique rows for cell arrays of strings/numbers.

n = size(A,1);
keys = cell(n,1);

for i = 1:n
    rowparts = cell(1,size(A,2));
    for j = 1:size(A,2)
        val = A{i,j};
        if isnumeric(val)
            rowparts{j} = num2str(val);
        else
            rowparts{j} = char(val);
        end
    end
    keys{i} = local_strjoin(rowparts, char(31));
end

[~, ia, ic] = unique(keys, 'stable');
C = A(ia,:);


function out = local_strjoin(parts, delimiter)
% Replacement for strjoin for older MATLAB versions.

if isempty(parts)
    out = '';
    return
end

out = parts{1};
for k = 2:numel(parts)
    out = [out, delimiter, parts{k}]; %#ok<AGROW>
end


function [mk, cc] = local_marker_color(i, Mark, colr)
% Select marker and color cyclically.

f1 = rem(i, length(Mark));
if f1 == 0
    f1 = length(Mark);
end

f2 = rem(i, size(colr,1));
if f2 == 0
    f2 = size(colr,1);
end

mk = Mark{f1};
cc = colr(f2,:);


function local_add_axes
% Add zero axes after data are plotted.

xx = axis;
plot(xx(1:2), [0 0], '--k', 'linewidth', 1.2)
plot([0 0], xx(3:4), '--k', 'linewidth', 1.2)
