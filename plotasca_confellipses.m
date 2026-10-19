function [hPts, hEll, out] = plotasca_confellipses(factormodel, sccomp, confLevel, whichScores, ascamodel)
% plotasca_confellipses Plot ASCA scores with confidence ellipses (per group)
%
% This function plots ASCA scores for a given design term and overlays
% confidence ellipses for each factor-level group in the selected 2D score space.
%
% The ellipse is a (1-confLevel) confidence region for the GROUP MEAN
% under multivariate normality with unknown covariance (Hotelling's T^2):
%
%   (z - mu)' * inv(S/n) * (z - mu) <= c
%   c = p*(n-1)/(n-p) * F_{p, n-p}(confLevel),  with p = 2
%
% INPUTS:
%   factormodel : field of ASCA model corresponding to the design term
%
% OPTIONAL INPUTS:
%   sccomp      : 1x2 vector of component indices to plot (default [1 2] if possible)
%   confLevel   : confidence level in (0,1), e.g. 0.95 (default 0.95)
%   whichScores : 'Scores' (default) or 'ScoreswithRes'
%   ascamodel   : full ASCA model (optional, to build legend labels)
%
% OUTPUTS:
%   hPts  : handles to point plots (one per group)
%   hEll  : handles to ellipse plots (one per group)
%   out   : struct with ellipse parameters (mu, S, n, c, coords)
%
% Author: (you can put your name)
% Inspired by Federico Marini's plotascascores grouping logic

% ---------------- Defaults ----------------
if nargin < 2 || isempty(sccomp)
    sccomp = 1:min(size(factormodel.SCA.Model.Scores,2), 2);
end
if numel(sccomp) < 2
    error('plotasca_confellipses requires TWO components in sccomp, e.g. [1 2].');
end
if nargin < 3 || isempty(confLevel)
    confLevel = 0.95;
end
if confLevel <= 0 || confLevel >= 1
    error('confLevel must be between 0 and 1 (e.g. 0.95).');
end
if nargin < 4 || isempty(whichScores)
    whichScores = 'Scores'; % or 'ScoreswithRes'
end

% Legend (optional, consistent with your plotascascores)
haveLegend = (nargin >= 5) && ~isempty(ascamodel) && isstruct(ascamodel);
legendLabels = {};
if haveLegend
    try
        legendLabels = local_build_legend_labels(factormodel, ascamodel);
        haveLegend = ~isempty(legendLabels);
    catch
        haveLegend = false;
    end
end

% ---------------- Extract score matrices ----------------
if ~isfield(factormodel.SCA.Model, 'Scores')
    error('factormodel.SCA.Model.Scores not found.');
end

haveProj = isfield(factormodel.SCA.Model, 'ScoreswithRes');

% Always keep mean scores available
Tmean = factormodel.SCA.Model.Scores;

% Decide which score set drives the ellipse/statistics
if strcmpi(whichScores,'ScoreswithRes')
    if ~haveProj
        error('ScoreswithRes requested but factormodel.SCA.Model.ScoreswithRes not found.');
    end
    Tell = factormodel.SCA.Model.ScoreswithRes;  % ellipse computed on this
else
    Tell = factormodel.SCA.Model.Scores;         % ellipse computed on this
end

if max(sccomp) > size(Tell,2)
    error('sccomp exceeds number of available components.');
end

T2_ell  = Tell(:,  sccomp(1:2));    % used for ellipses
T2_mean = Tmean(:, sccomp(1:2));    % used for filled markers

if strcmpi(whichScores,'ScoreswithRes')
    T2_res = factormodel.SCA.Model.ScoreswithRes(:, sccomp(1:2)); % empty markers
else
    T2_res = [];
end

% Grouping index from design matrix rows
[~,~,lindex] = unique(factormodel.DesignMatrix,'rows','stable');
ng = max(lindex);

% Style settings
Mark = {'o','s','d','^','v','>','<','p','h','+','*','.','x'};
colr = lines(7);

figure; hold on;
hPts = gobjects(ng,1);          % handle for filled mean points
hPtsRes = gobjects(ng,1);       % handle for empty projected points (if applicable)
hEll = gobjects(ng,1);

out = struct();
out.confLevel = confLevel;
out.whichScores = whichScores;
out.sccomp = sccomp(1:2);
out.groups = cell(ng,1);

p = 2;
theta = linspace(0, 2*pi, 200);
unitCircle = [cos(theta); sin(theta)];

for g = 1:ng
    idx = (lindex == g);

    % group styles
    f1 = rem(g,13); if f1==0, f1=13; end
    f2 = rem(g,7);  if f2==0, f2=7;  end
    thisColor = colr(f2,:);
    thisMark  = Mark{f1};

    % ---- Filled markers: mean scores (always)
    Ym = T2_mean(idx,:);
    hPts(g) = plot(Ym(:,1), Ym(:,2), thisMark, ...
        'MarkerFaceColor', thisColor, ...
        'MarkerEdgeColor', thisColor, ...
        'MarkerSize', 10);

    % ---- Empty markers: ScoreswithRes (only if requested)
    if strcmpi(whichScores,'ScoreswithRes')
        Yr = T2_res(idx,:);
        hPtsRes(g) = plot(Yr(:,1), Yr(:,2), thisMark, ...
            'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', thisColor, ...
            'MarkerSize', 10);
    end

    % ---- Ellipse computed from selected score set (Tell)
    Y = T2_ell(idx,:);
    n = size(Y,1);

    if n <= p
        out.groups{g} = struct('n',n,'mu',mean(Y,1), ...
            'S',nan(2),'c',nan,'ellipse',nan(2,numel(theta)));
        continue;
    end

    mu = mean(Y, 1);
    S  = cov(Y, 1);              % normalized by n
    c  = p*(n-1)/(n-p) * finv(confLevel, p, n-p);

    Smean = (S / n);
    Smean = (Smean + Smean')/2;

    [R,flag] = chol(Smean, 'lower');
    if flag ~= 0
        [V,D] = eig(Smean);
        D = max(D,0);
        R = V*sqrt(D);
    end

    ell = mu(:) + sqrt(c) * (R * unitCircle);

    % IMPORTANT: set ellipse color to match group points
    hEll(g) = plot(ell(1,:), ell(2,:), '-', 'LineWidth', 1.8, 'Color', thisColor);

    out.groups{g} = struct('n',n,'mu',mu,'S',S,'c',c,'ellipse',ell);
end

% Axes cosmetics
xx = axis;
plot(xx(1:2), [0 0], '--k', 'linewidth', 1.2);
plot([0 0], xx(3:4), '--k', 'linewidth', 1.2);

xlabel(['Scores on SC', num2str(sccomp(1))], 'fontsize', 12, 'fontweight', 'bold');
ylabel(['Scores on SC', num2str(sccomp(2))], 'fontsize', 12, 'fontweight', 'bold');
set(gca, 'linewidth', 1.5, 'fontsize', 12, 'fontweight', 'bold');

if haveLegend
    legend([hPts(1), hPtsRes(1), hEll(1)], {'Scores (mean)','ScoreswithRes (proj)','Conf. ellipse'}, 'Location','best');
end

title(sprintf('ASCA %s with %.1f%% confidence ellipses (group mean)', whichScores, 100*confLevel), ...
      'fontsize', 12, 'fontweight', 'bold');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function labels = local_build_legend_labels(factormodel, model)
% Same idea as in plotascascores.m: builds labels from original DesignInput.

labels = {};

if isfield(model,'DesignInput')
    D0 = model.DesignInput;
elseif isfield(model,'Design')
    D0 = model.Design;
else
    return;
end

terms = {};
if isfield(factormodel,'Members')
    terms = factormodel.Members;
elseif isfield(factormodel,'DesignTerms')
    dt = factormodel.DesignTerms;
    if isnumeric(dt)
        terms = arrayfun(@(x) char(64+x), dt(:)', 'UniformOutput', false);
    end
end
if isempty(terms)
    if isfield(factormodel,'EffectLabel')
        terms = {factormodel.EffectLabel};
    else
        return;
    end
end

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
fidx = double(letters) - 64;

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
end

function s = local_level_string(D0, rowIdx, fidx, model)
parts = cell(1, numel(fidx));
for k = 1:numel(fidx)
    j = fidx(k);
    fname = char(64+j);

    if isfield(model,'DesignInfo') && isfield(model.DesignInfo,'Factors') && numel(model.DesignInfo.Factors) >= j
        if isfield(model.DesignInfo.Factors(j), 'Name') && ~isempty(model.DesignInfo.Factors(j).Name)
            fname = model.DesignInfo.Factors(j).Name;
        end
    end

    if isnumeric(D0)
        lv = num2str(D0(rowIdx, j));
    elseif iscell(D0)
        lv = D0{rowIdx, j};
        if isnumeric(lv), lv = num2str(lv); else, lv = char(string(lv)); end
    else
        lv = char(string(D0(rowIdx, j)));
    end

    parts{k} = sprintf('%s=%s', fname, lv);
end
s = strjoin(parts, ', ');
end
