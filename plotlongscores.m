
function plotlongscores(varargin)
%PLOTLONGSCORES Longitudinal scores plots for ASCA / ASCA-flex models.
%
% Usage (recommended):
%   plotlongscores(fullModel, termLabel, 'timefactor', i, 'groupfactor', j, 'comp', c)
%   plotlongscores(termModel, 'model', fullModel, 'timefactor', i, 'groupfactor', j, 'comp', c)
%
% termLabel can be an internal label (e.g. 'A','B_AB') or display label (e.g. 'B+AB').
% timefactor and groupfactor refer to columns of fullModel.DesignInput/Design (original factor columns).
%
% If the time-factor values are strings that look numeric (e.g. '40','50'), this function
% will convert them to numeric time points for tscale. Otherwise it treats time as categorical.
%
% Backward compatible:
%   plotlongscores(fullModel, termLabel)
%   plotlongscores(fullModel, termLabel, tscale)
%   plotlongscores(fullModel, termLabel, tscale, flab)

% -----------------------------
% Legacy positional parsing
% -----------------------------
a1 = varargin{1};
ascamodel = [];
termLabel = '';
tscale_in = [];
flab_in = [];

if isstruct(a1) && (isfield(a1,'Blocks') || isfield(a1,'Design') || isfield(a1,'DesignInput'))
    % full model first
    ascamodel = a1;
    if nargin >= 2 && (ischar(varargin{2}) || isstring(varargin{2}))
        termLabel = char(varargin{2});
    end
    rest = varargin(3:end);
else
    % term model first
    termModel = a1;
    rest = varargin(2:end);
end

% -----------------------------
% Name-value parsing
% -----------------------------
timefactor = [];
groupfactor = [];
comp = 1;
useRes = false;
fullModelNV = [];

k = 1;
while k <= numel(rest)
    if ischar(rest{k}) || isstring(rest{k})
        key = lower(char(rest{k}));
        if strcmp(key,'timefactor')
            timefactor = rest{k+1}; k = k + 2; continue;
        elseif strcmp(key,'groupfactor')
            groupfactor = rest{k+1}; k = k + 2; continue;
        elseif strcmp(key,'comp') || strcmp(key,'component')
            comp = rest{k+1}; k = k + 2; continue;
        elseif strcmp(key,'usewithres') || strcmp(key,'useres')
            useRes = logical(rest{k+1}); k = k + 2; continue;
        elseif strcmp(key,'model')
            fullModelNV = rest{k+1}; k = k + 2; continue;
        else
            k = k + 1; continue;
        end
    else
        k = k + 1;
    end
end

% Resolve full model / term model
if isempty(ascamodel)
    if ~isempty(fullModelNV)
        ascamodel = fullModelNV;
    else
        error('When first input is a term model, pass the full model via ''model'', e.g. plotlongscores(termModel,''model'',fullModel,...)');
    end
end

% Resolve term model
if exist('termModel','var') && ~isempty(termModel)
    % already provided
    internalLabel = '';
else
    if isempty(termLabel)
        error('termLabel is required when first input is the full model.');
    end
    [internalLabel, termModel] = local_get_term_model(ascamodel, termLabel);
end

if isempty(internalLabel)
    % Try to infer internal label from termModel if possible (rare)
    internalLabel = local_guess_internal_label(ascamodel, termModel);
end

% Determine default timefactor/groupfactor if not provided:
% If the term has DesignTerms, pick first as timefactor
if isempty(timefactor)
    if isfield(termModel,'DesignTerms') && ~isempty(termModel.DesignTerms)
        timefactor = termModel.DesignTerms(1);
    else
        timefactor = 1;
    end
end

% Scores
if useRes && isfield(termModel,'SCA') && isfield(termModel.SCA,'Model') && isfield(termModel.SCA.Model,'ScoreswithRes')
    T = termModel.SCA.Model.ScoreswithRes;
else
    T = termModel.SCA.Model.Scores;
end
if comp < 1 || comp > size(T,2)
    error('Requested comp=%d but term has %d components.', comp, size(T,2));
end
score = T(:, comp);

% Factor columns (per-sample) from model
timeCol = local_get_factor_column(ascamodel, timefactor);
[groupCol, haveGroup] = local_get_group_column(ascamodel, groupfactor);

% Build time scale + labels (numeric strings -> numeric)
[tlevels, tpos, tticks, tticklab] = local_build_time_scale(timeCol, tscale_in);

% Build group levels + legend labels
if haveGroup
    [glevels, glab] = local_build_group_levels(groupCol, flab_in, ascamodel, groupfactor);
else
    glevels = {''};
    glab = {''};
end

% Aggregate mean score per (group,time)
Y = nan(numel(glevels), numel(tlevels));
Ys = nan(numel(glevels), numel(tlevels));

for gi = 1:numel(glevels)
    for ti = 1:numel(tlevels)
        idx = local_match_level(timeCol, tlevels{ti});
        if haveGroup
            idx = idx & local_match_level(groupCol, glevels{gi});
        end
        if any(idx)
            Y(gi,ti) = mean(score(idx));
            if numel(idx) > 1
                Ys(gi,ti) = std(score(idx)) ./ sqrt(numel(idx));
            else
                Ys(gi,ti) = 0;
            end
        end
    end
end

% Plot
figure;
hold on;
h = gobjects(numel(glevels),1);
for gi = 1:numel(glevels)
    h(gi) = errorbar(tpos, Y(gi,:), Ys(gi,:), '-o', 'LineWidth', 1.5);
end
hold off;
grid on;
xlabel(local_factor_name(ascamodel, timefactor));
ylabel(sprintf('Score comp %d', comp));
set(gca, 'fontsize', 15, 'fontweight', 'bold', 'linewidth',2)

title(local_term_title(ascamodel, internalLabel));

set(gca,'XTick', tticks);
set(gca,'XTickLabel', tticklab);

if haveGroup
    legend(h, glab, 'Interpreter','none', 'Location','best');
end

end

% -----------------------------
% Helpers
% -----------------------------
function [internalLabel, termModel] = local_get_term_model(m, termLabel)
% Resolve internal/display label to model field.
q = strrep(strtrim(termLabel), ' ', '');
q = strrep(q, '+', '_');

% direct
cand = ['X' q];
if isfield(m, cand)
    internalLabel = q;
    termModel = m.(cand);
    return;
end

% if Blocks exist, search
if isfield(m,'Blocks') && ~isempty(m.Blocks)
    for k = 1:numel(m.Blocks)
        il = strrep(m.Blocks(k).InternalLabel,' ','');
        dl = strrep(m.Blocks(k).DisplayLabel,' ','');
        dl = strrep(dl,'+','_');
        if strcmpi(q, il) || strcmpi(q, dl)
            internalLabel = m.Blocks(k).InternalLabel;
            termModel = m.(['X' internalLabel]);
            return;
        end
    end
end

% last ditch: any X* field
fn = fieldnames(m);
for k = 1:numel(fn)
    if startsWith(fn{k},'X') && strcmpi(fn{k}(2:end), q)
        internalLabel = fn{k}(2:end);
        termModel = m.(fn{k});
        return;
    end
end

error('Could not resolve term label "%s" to a model field.', termLabel);
end

function il = local_guess_internal_label(m, termModel)
il = '';
if ~isstruct(termModel); return; end
fn = fieldnames(m);
xfields = fn(startsWith(fn,'X'));
for k=1:numel(xfields)
    if isequaln(m.(xfields{k}), termModel)
        il = xfields{k}(2:end);
        return;
    end
end
end

function col = local_get_factor_column(m, fidx)
% Prefer DesignInput; accept DesignInfo cell matrix too.
if isfield(m,'DesignInput') && ~isempty(m.DesignInput)
    D = m.DesignInput;
elseif isfield(m,'DesignInfo') && iscell(m.DesignInfo) && isfield(m,'Design') && size(m.DesignInfo,1) == size(m.Design,1)
    D = m.DesignInfo; % per-sample cell matrix
elseif isfield(m,'Design') && ~isempty(m.Design)
    D = m.Design;
else
    error('Model does not contain DesignInput/Design/DesignInfo needed for plotting.');
end
if fidx < 1 || fidx > size(D,2)
    error('Factor index %d out of range (design has %d columns).', fidx, size(D,2));
end
col = D(:, fidx);
end

function [gcol, haveGroup] = local_get_group_column(m, gidx)
if isempty(gidx)
    haveGroup = false;
    gcol = [];
    return;
end
haveGroup = true;
gcol = local_get_factor_column(m, gidx);
end

function [tlevels, tpos, tticks, tticklab] = local_build_time_scale(timeCol, tscale_in)
% timeCol: per-sample column values (numeric / cellstr / string)
% Returns:
%   tlevels: cell array of unique level tokens
%   tpos: numeric x positions for each level
%   tticks/tticklab: axis ticks & labels

% If user provided tscale explicitly (legacy), treat time as positions 1..k unless numeric vector same length.
if ~isempty(tscale_in)
    if isnumeric(tscale_in)
        tpos = tscale_in(:)';
        tticks = tpos;
        tticklab = arrayfun(@num2str, tpos, 'UniformOutput', false);
        % levels taken from unique stable in data
    else
        % categorical labels
        tticklab = tscale_in;
        tpos = 1:numel(tticklab);
        tticks = tpos;
    end
else
    tpos = [];
    tticks = [];
    tticklab = {};
end

% Build unique levels in stable order from timeCol
[tlevels_raw, ia] = local_unique_stable(timeCol);
tlevels = tlevels_raw; % tokens as stored

% Try numeric conversion if timeCol is cell/string AND all unique levels parse
if isempty(tscale_in)
    nums = local_try_num(tlevels_raw);
    if ~isempty(nums)
        % Use numeric times; sort ascending for longitudinal x-axis
        [numsSorted, ord] = sort(nums);
        tlevels = tlevels_raw(ord);
        tpos = numsSorted(:)';
        tticks = tpos;
        tticklab = arrayfun(@num2str, numsSorted(:), 'UniformOutput', false);
    else
        % categorical
        tpos = 1:numel(tlevels_raw);
        tticks = tpos;
        tticklab = local_to_cellstr(tlevels_raw);
    end
end
end

function nums = local_try_num(levels)
% levels: cell tokens; return numeric if all parse and none NaN
nums = [];
try
    s = local_to_cellstr(levels);
    v = str2double(s);
    if all(~isnan(v))
        nums = v(:);
    end
catch
    nums = [];
end
end

function [levels, labels] = local_build_group_levels(groupCol, flab_in, model, groupfactor)
% Returns group level tokens and legend labels.
[levels, ~] = local_unique_stable(groupCol);

if ~isempty(flab_in)
    labels = flab_in;
    if isstring(labels), labels = cellstr(labels); end
    return;
end

% If model has DesignInfo.Factors with OriginalLevels, use them to label
labels = local_to_cellstr(levels);

% If original levels are numeric strings, keep as-is in legend.
% Optionally prefix factor name
fname = local_factor_name(model, groupfactor);
for i=1:numel(labels)
    if ~isempty(fname)
        labels{i} = [fname '=' labels{i}];
    end
end
end

function idx = local_match_level(col, levelToken)
% Match per-sample column entries to a given token (handles numeric/cell/string)
if isnumeric(col)
    idx = (col == levelToken);
elseif iscell(col)
    idx = ismember(col, levelToken);
elseif isstring(col)
    idx = (col == string(levelToken));
else
    idx = false(size(col,1),1);
end
end

function [u, ia] = local_unique_stable(col)
% Returns unique levels in stable order, as cell tokens.
if isnumeric(col)
    [uu, ia] = unique(col, 'stable');
    u = num2cell(uu);
elseif isstring(col)
    colc = cellstr(col);
    [u, ia] = unique(colc, 'stable');
    u = u(:);
elseif iscell(col)
    [u, ia] = unique(col, 'stable');
    u = u(:);
else
    error('Unsupported design column type for plotting.');
end
end

function c = local_to_cellstr(levels)
if iscell(levels)
    % could contain numeric
    if ~isempty(levels) && isnumeric(levels{1})
        c = cellfun(@num2str, levels, 'UniformOutput', false);
    else
        c = levels;
    end
elseif isstring(levels)
    c = cellstr(levels);
else
    c = cellfun(@num2str, num2cell(levels), 'UniformOutput', false);
end
end

function s = local_factor_name(m, fidx)
s = sprintf('F%d', fidx);
if isfield(m,'DesignInfo')
    try
        if isstruct(m.DesignInfo) && isfield(m.DesignInfo,'Factors') && numel(m.DesignInfo.Factors) >= fidx
            if isfield(m.DesignInfo.Factors(fidx),'Name') && ~isempty(m.DesignInfo.Factors(fidx).Name)
                s = m.DesignInfo.Factors(fidx).Name;
            end
        end
    catch
        % ignore
    end
end
end

function t = local_term_title(m, internalLabel)
t = internalLabel;
if isfield(m,'Blocks') && ~isempty(m.Blocks)
    for k=1:numel(m.Blocks)
        if strcmp(m.Blocks(k).InternalLabel, internalLabel)
            t = m.Blocks(k).DisplayLabel;
            return;
        end
    end
end
end
