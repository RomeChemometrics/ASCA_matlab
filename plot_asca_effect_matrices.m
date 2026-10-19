function out = plot_asca_effect_matrices(model, varargin)
% plot_asca_effect_matrices  Plot all ASCA effect matrices in subplots
%
% out = plot_asca_effect_matrices(model, 'Name',Value,...)
%
% This function expects a model struct from asca_flex.m where effects are in
% fields like model.XA, model.XB, model.XAB, model.XRes, etc.
% Each effect field contains .EffectMatrix and .EffectLabel (usually).
%
% -------------------- Options --------------------
% 'Scale'     : scaling applied to columns before plotting (default 'none')
%              - 'none'   : no scaling
%              - 'zscore' : mean-center and divide by std (per column)
%              - 'pareto' : mean-center and divide by sqrt(std)
%              - numeric vector (1xK or Kx1): multiply each column by this factor
% 'XLabel'    : string for x-axis label in each subplot (default 'Variables')
% 'VarNames'  : cellstr of length K for x tick labels (default [])
% 'ShowTicks' : true/false show x ticks (default: true if VarNames provided, else false)
% 'CLim'      : [] or [cmin cmax] for color scaling (default [])
% 'Colormap'  : colormap name or Mx3 array (default: parula)
% 'MaxCols'   : max number of x tick labels to show (default 20)
%
% -------------------- Output --------------------
% out.Effects: struct array with effect metadata and plotted matrix
%
% Example:
%   plot_asca_effect_matrices(model, 'Scale','zscore', 'XLabel','m/z');
%   plot_asca_effect_matrices(model, 'Scale',ones(1,K), 'VarNames',mzLabels);

% ---------- Parse inputs ----------
p = inputParser;
p.addParameter('Scale', 'none');
p.addParameter('XLabel', 'Variables', @(s) ischar(s) || isstring(s));
p.addParameter('VarNames', {}, @(x) iscell(x) || isstring(x) || isempty(x));
p.addParameter('ShowTicks', [], @(x) islogical(x) || isempty(x));
p.addParameter('CLim', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
p.addParameter('Colormap', parula, @(x) isnumeric(x) || ischar(x) || isstring(x));
p.addParameter('MaxCols', 20, @(x) isnumeric(x) && isscalar(x) && x>=2);
p.addParameter('PlotType','image',@(s)ischar(s)||isstring(s));
p.addParameter('XValues', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.parse(varargin{:});
opt = p.Results;

% normalize var names
if isstring(opt.VarNames), opt.VarNames = cellstr(opt.VarNames); end

% ---------- Find effect fields ----------
fn = fieldnames(model);
effects = struct('Field',{},'Label',{},'Matrix',{});

for i = 1:numel(fn)
    f = fn{i};
    if ~startsWith(f,'X');  continue; end
    if ~isstruct(model.(f)); continue; end
    if ~isfield(model.(f),'EffectMatrix'); continue; end

    M = model.(f).EffectMatrix;
    if ~isnumeric(M) || ndims(M) ~= 2
        continue
    end

    if isfield(model.(f),'EffectLabel') && ~isempty(model.(f).EffectLabel)
        lab = model.(f).EffectLabel;
    else
        lab = f;
    end

    effects(end+1).Field = f; %#ok<AGROW>
    effects(end).Label  = string(strtrim(lab));
    effects(end).Matrix = M;
end

if isempty(effects)
    error('No fields with .EffectMatrix found. Are you passing the output of asca_flex()?');
end

% Sort: put Mean first if present, then main effects, then interactions, then Res
labels = lower(string({effects.Label}));

isMean = contains(labels,"mean");
isRes  = contains(labels,"res");

% put mean first, residual last, keep the rest in original order
idxSort = [ find(isMean) , ...
    find(~isMean & ~isRes) , ...
    find(isRes) ];

effects = effects(idxSort);

% ---------- Determine subplot grid ----------
nE = numel(effects);
nR = floor(sqrt(nE));
nC = ceil(nE/nR);
if nR*nC < nE
    nR = ceil(nE/nC);
end

% ---------- Plot ----------
figure;
set(gcf,'Color','w','units', 'normalized', 'position',[0.1 0.1 0.8 0.8]);

out = struct();
out.Effects = struct([]);

for e = 1:nE
    M = effects(e).Matrix;  % N x K
    [N,K] = size(M);

    % Validate VarNames length (optional)
    if ~isempty(opt.VarNames) && numel(opt.VarNames) ~= K
        error('VarNames length (%d) must match number of columns K=%d.', numel(opt.VarNames), K);
    end


    % Scaling
    Mplot = apply_scaling(M, opt.Scale);

    % X-axis
    xv = opt.XValues;
    if ~isempty(xv)
        if numel(xv) ~= K
            error('XValues must have the same length as the number of variables (%d).', K);
        end
        xv = xv(:)';   % row
    end

    % Subplot
    ax = subplot(nR, nC, e);

    plotType = lower(string(opt.PlotType));

    switch plotType

        case "image"

            if isempty(xv)
                imagesc(Mplot);
            else
                imagesc(xv, 1:size(Mplot,1), Mplot);
            end
            axis tight;

        case "line"

            if isempty(xv)
                plot(Mplot','LineWidth',2);
            else
                plot(xv, Mplot','LineWidth',2);
            end
            axis tight;
            yline(0,'k:');

        otherwise
            error('PlotType must be ''image'' or ''line''.');
    end

    title(effects(e).Label, 'FontSize', 15, 'FontWeight', 'bold');

    % Axis labels
    xlabel(opt.XLabel, 'FontSize', 15, 'FontWeight', 'bold');
    ylabel('Samples',  'FontSize', 15, 'FontWeight', 'bold');

    % Colormap + limits
    % Colormap (robust across MATLAB releases)
    cm = opt.Colormap;

    if isstring(cm)
        cm = char(cm);   % convert "parula" -> 'parula'
    end

    if (ischar(cm) || isstring(opt.Colormap))
        % colormap name
        try
            colormap(cm);         % works broadly
        catch
            colormap(gca, cm);    % fallback
        end
    elseif isnumeric(cm) && size(cm,2) == 3
        % numeric Mx3 colormap
        try
            colormap(cm);
        catch
            colormap(gca, cm);
        end
    else
        error('Colormap must be a colormap name (char/string) or an Mx3 numeric array.');
    end
    if ~isempty(opt.CLim)
        clim(ax, opt.CLim);
    end

    % ticks
    if isempty(opt.ShowTicks)
        showTicks = ~isempty(opt.VarNames);
    else
        showTicks = opt.ShowTicks;
    end

    if showTicks

        % ---- case 1: variable names available (categorical axis)
        if isempty(xv) && ~isempty(opt.VarNames)

            maxCols = opt.MaxCols;
            if K <= maxCols
                xticks(1:K);
                xticklabels(opt.VarNames);
            else
                step = ceil(K/maxCols);
                idx  = 1:step:K;
                xticks(idx);
                xticklabels(opt.VarNames(idx));
            end
            xtickangle(45);

            % ---- case 2: numeric x-axis or no labels provided
        else
            % just let MATLAB choose reasonable numeric ticks
            xticks('auto');
        end

    else
        xticks([]);
    end

    % Styling requested
    set(ax, 'FontSize', 15, 'FontWeight', 'bold', 'LineWidth', 2);

    % Colorbar per subplot (keeps it simple)
    if plotType == "image"
        colorbar;
    end
    out.Effects(e).Field = effects(e).Field;
    out.Effects(e).Label = effects(e).Label;
    out.Effects(e).N     = N;
    out.Effects(e).K     = K;
    out.Effects(e).Scale = opt.Scale;
end

end

% =================== Helper: scaling ===================
function Mout = apply_scaling(M, scaleOpt)

if isnumeric(scaleOpt)
    s = scaleOpt(:)'; % 1xK
    if numel(s) ~= size(M,2)
        error('Numeric Scale must have length K=%d.', size(M,2));
    end
    Mout = M .* repmat(s, size(M,1), 1);
    return
end

scaleOpt = lower(string(scaleOpt));

switch scaleOpt
    case "none"
        Mout = M;

    case "zscore"
        mu = mean(M,1,'omitnan');
        sd = std(M,0,1,'omitnan');
        sd(sd==0) = 1;
        Mout = (M - mu) ./ sd;

    case "pareto"
        mu = mean(M,1,'omitnan');
        sd = std(M,0,1,'omitnan');
        sdp = sqrt(sd);
        sdp(sdp==0) = 1;
        Mout = (M - mu) ./ sdp;

    otherwise
        error('Unknown Scale option: %s. Use none|zscore|pareto or a numeric vector.', scaleOpt);
end

end
