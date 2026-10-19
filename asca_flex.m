function model=asca_flex(varargin)
% asca Calculates an ASCA model on designed data
%  The function asca calculates an ANOVA-Simultaneous Component Analysis (ASCA)
%  model for the analysis of multivariate data resulting from a designed
%  experiment.
%
%  INPUT:
%            X = Matrix collecting the experimental data
%    desmatrix = Matrix of the experimental design (levels should be
%                coded numerically (even if not necessarily
%                consequentially).
%  OPTIONAL INPUT
%          opt = Structure collecting the options for training the model
%                Contains the following fields.
%  reducedmodel: [{'standard'}] Defines the calculation of the
%                reduced ANOVA model for modeling and testing. Alternatively,
%                it can be a cell array, each component containing a vector of
%                terms to be removed to calculate customly the residuals under
%                the reduced model
%      permtest: [{'on'} | 'off'] Decides if permutation tests have to be
%                performed or not for factor significance
%         nperm: Number of permutations (default = 10000)
%         plots: [{'on'} | 'off'] Regulates plotting
%     bootstrap: [{'off'} | 'all' | 'signif'] Decides if bootstrap of SCA loadings
%                has to be performed or not for evaluation of variables'
%                contributions. If 'signif' is selected, bootstrap is only
%                applied to the SCA model of the terms identified as
%                significant by the permutation test.
%         nboot: Number of bootstrap repetitions (default = 10000)
%      bootsave: ['all' | {'confint'} | 'signvars'] Regulates the amount of
%                bootstrap results saved in the model.
%       preproc: defines the pretreatment of X matrix before ASCA analysis
%                (default is {'none'})
%         confl: Confidence level for bootstrap confidence interval
%
%
%       Default options can be obtained by typing:
%       opt=asca('options');
%
% . OUTPUT:
%        model = Structure with the final ASCA model
%
%   I/O:
%      model=asca(X, desmatrix, opt);
%      model=asca(X, desmatrix);
%
%
% Author: Federico Marini
% Version: 14/05/2020



if nargin==1 && strcmp(varargin{1}, 'options')

    opt.preproc='none';
    opt.reducedmodel='standard';
    opt.permtest='on';
    opt.permfacts='all';
    opt.permmethod = 'freedman-lane';
    opt.nperm=10000;
    opt.sca_source = 'typeiii';  %alternative is 'typei'
    opt.bootstrap='off';
    opt.bootmatrix='original';
    opt.bootsave='confint';
    opt.nboot=1000;
    opt.plots='on';
    opt.confl=0.95;
    opt.maxorder=inf;           % maximum interaction order to include (inf = full)
    opt.includeterms='all';      % 'all' or cell array of term labels, e.g. {'A','B','AB'}
    opt.termgroups={};           % optional grouping, e.g. {{'C','AC'}} or {'C+AC'}

    model=opt;
    return

elseif nargin==2
    X=varargin{1};
    desmatrix=varargin{2};
    opt.preproc='none';
    opt.reducedmodel='standard';
    opt.permtest='on';
    opt.permfacts='all';
    opt.permmethod = 'freedman-lane';
    opt.nperm=10000;
    opt.sca_source = 'typeiii';  %alternative is 'typei'
    opt.bootstrap='off';
    opt.bootmatrix='original';
    opt.bootsave='confint';
    opt.nboot=1000;
    opt.plots='on';
    opt.confl=0.95;
    opt.maxorder=inf;           % maximum interaction order to include (inf = full)
    opt.includeterms='all';      % 'all' or cell array of term labels, e.g. {'A','B','AB'}
    opt.termgroups={};           % optional grouping, e.g. {{'C','AC'}} or {'C+AC'}

elseif nargin==3
    X=varargin{1};
    desmatrix=varargin{2};
    opt=varargin{3};
end


%preprocessing
ntot=size(X,1);

switch opt.preproc
    case 'none'
        Xp=X;
    case 'auto'
        Xp=X./repmat(std(X), ntot, 1);
    case 'pareto'
        Xp=X./repmat(sqrt(std(X)),ntot, 1);
end

model.Xdata.PreprData=Xp;
model.Xdata.TotalSSQ=sum(sum(Xp.^2));

% --- Normalize design input (allow non-consecutive numeric levels or labels) ---
[desmatrix_coded, designInfo] = normalize_design_input(desmatrix);
model.DesignInput = desmatrix;          % original user input (numeric or cell labels)
model.DesignInfo  = designInfo;         % mapping per factor
model.Design      = desmatrix_coded;    % numeric coded (1..k per factor)



dmain=gendmat(desmatrix_coded);
[dmatrices,desterms,deslabels, desorder]=createdesign(dmain);


% ------------------------------------------------------------
% Term selection / truncation / grouping (optional)
% ------------------------------------------------------------
if ~isfield(opt,'maxorder');    opt.maxorder = inf; end
if ~isfield(opt,'includeterms'); opt.includeterms = 'all'; end
if ~isfield(opt,'termgroups');  opt.termgroups = {}; end

% Work with trimmed labels
tlabels = cellfun(@strtrim, deslabels, 'UniformOutput', false);

% 1) Truncate by interaction order (Mean is always kept)
keep = true(length(tlabels),1);
keep(1) = true;
if ~isinf(opt.maxorder)
    keep = keep & (desorder <= opt.maxorder);
    keep(1) = true;
end

% 2) Include only specified terms (optional)
if ~(ischar(opt.includeterms) && strcmpi(opt.includeterms,'all')) && ~isempty(opt.includeterms)
    inc = opt.includeterms;
    if ischar(inc); inc = {inc}; end
    inc = cellfun(@canonical_label, inc, 'UniformOutput', false);
    keep = keep & ismember(tlabels, inc);
    keep(1) = true; % always include Mean
end

% Apply selection
dmatrices = dmatrices(keep);
desterms  = desterms(keep);
desorder  = desorder(keep);
deslabels = deslabels(keep);
tlabels   = tlabels(keep);

% 3) Group terms together (optional)
groupInfo = struct('Label',{},'DisplayLabel',{},'Members',{},'MemberIdx',{});
if ~isempty(opt.termgroups)
    used = false(length(tlabels),1); % track terms that get grouped
    new_dmat = {};
    new_terms = {};
    new_labels = {};
    new_order = [];

    for g = 1:length(opt.termgroups)
        grp = opt.termgroups{g};
        if ischar(grp)
            grp = strsplit(grp, '+');
        end
        if isstring(grp)
            grp = cellstr(grp);
        end

        grp = cellfun(@strtrim, grp, 'UniformOutput', false);
        grp = cellfun(@canonical_label, grp, 'UniformOutput', false);

        idx = find(ismember(tlabels, grp));
        if isempty(idx)
            continue
        end

        used(idx) = true;

        % concatenate design matrices of grouped terms
        Dg = [];
        for kk = 1:length(idx)
            Dg = [Dg, dmatrices{idx(kk)}];
        end

        % grouped label as provided order (canonicalized)
        glab_display = strjoin(grp, '+');
        glab = strjoin(grp, '_'); % internal label for struct field names

        new_dmat{end+1}   = Dg;
        new_terms{end+1}  = unique([desterms{idx}]);
        new_labels{end+1} = glab;
        new_order(end+1)  = max(desorder(idx));

        groupInfo(end+1).Label     = glab;
        groupInfo(end).DisplayLabel = glab_display;
        groupInfo(end).Members     = grp;
        groupInfo(end).MemberIdx   = idx;
    end

    % keep ungrouped terms
    keep2 = ~used;
    dmatrices = dmatrices(keep2);
    desterms  = desterms(keep2);
    desorder  = desorder(keep2);
    deslabels = deslabels(keep2);
    tlabels   = tlabels(keep2);

    % append grouped terms
    if ~isempty(new_dmat)
        dmatrices = [dmatrices, new_dmat];
        desterms  = [desterms, new_terms];
        desorder  = [desorder; new_order(:)];
        deslabels = [deslabels, new_labels];
    end

    % Re-sort (Mean first, then by order then label)
    tlabels = cellfun(@strtrim, deslabels, 'UniformOutput', false);
    if ~isempty(tlabels)
        isMean = strcmpi(tlabels,'Mean');
        idxMean = find(isMean,1);
        others = find(~isMean);
        if ~isempty(others)
            key = [desorder(others), double(char(tlabels(others)))];
            [~, ord] = sortrows(key);
            others = others(ord);
        end
        if ~isempty(idxMean)
            newidx = [idxMean; others(:)];
        else
            newidx = others(:);
        end
        dmatrices = dmatrices(newidx);
        desterms  = desterms(newidx);
        desorder  = desorder(newidx);
        deslabels = deslabels(newidx);
        tlabels   = tlabels(newidx);
    end
end

model.TermSelection.maxorder    = opt.maxorder;
model.TermSelection.includeterms = opt.includeterms;
model.TermSelection.termgroups  = opt.termgroups;
model.TermSelection.groups      = groupInfo;
% Human-readable labels (for grouped terms we prefer "A+AB" style)
model.TermDisplayLabels = deslabels;
if ~isempty(groupInfo)
    for gg = 1:length(groupInfo)
        idxg = find(strcmp(model.TermDisplayLabels, groupInfo(gg).Label));
        if ~isempty(idxg)
            model.TermDisplayLabels{idxg} = groupInfo(gg).DisplayLabel;
        end
    end
end


% refresh labels used downstream
deslabels = tlabels;
model.TermLabels = deslabels;


% ------------------------------------------------------------
% Block metadata (post-selection / post-grouping)
% InternalLabel is used for struct field names; DisplayLabel is for reporting
% Members lists the original (canonical) terms included in the block.
% ------------------------------------------------------------
nBlocks = length(dmatrices);
Blocks = struct('InternalLabel',{},'DisplayLabel',{},'Members',{},'Order',{},'DesignMatrix',{},'ColIdx',{});
for bi = 1:nBlocks
    ilab = strtrim(model.TermLabels{bi});
    dlab = strtrim(model.TermDisplayLabels{bi});
    mem  = {ilab};
    if ~isempty(groupInfo)
        gi = find(strcmp({groupInfo.Label}, ilab), 1);
        if ~isempty(gi)
            mem = groupInfo(gi).Members;
        end
    end
    Blocks(bi).InternalLabel = ilab;
    Blocks(bi).DisplayLabel  = dlab;
    Blocks(bi).Members       = mem;
    Blocks(bi).Order         = desorder(bi);
    Blocks(bi).DesignMatrix  = dmatrices{bi};
    Blocks(bi).ColIdx        = [];
end
model.Blocks = Blocks;




Xd=Xp;

for i=1:length(dmatrices)
    Xeff=dmatrices{i}*pinv(dmatrices{i})*Xd;
    ssqEff=sum(sum(Xeff.^2));
    Xd=Xd-Xeff;
    l=strtrim(deslabels{i});
    eval(['model.X', l,'.EffectMatrix=Xeff;'])
    eval(['model.X', l,'.EffectSSQ=ssqEff;'])

    if i==1
        ssqtot=sum(sum(Xd.^2));
        model.Xdata.CenteredData=Xd;
        model.Xdata.CenteredSSQ=ssqtot;
    else
        expVar=100*(ssqEff/ssqtot);
        eval(['model.X', l,'.EffectExplVar=expVar;'])
    end

    eval(['model.X', l,'.DesignMatrix=dmatrices{i};'])
    eval(['model.X', l,'.DesignTerms=desterms{i};'])
    eval(['model.X', l,'.EffectLabel=strtrim(deslabels{i});'])
    eval(['model.X', l,'.TermOrder=desorder(i);'])

end


model.XRes.EffectMatrix=Xd;
model.XRes.EffectSSQ=sum(sum(Xd.^2));
model.XRes.EffectExplVar=100*(model.XRes.EffectSSQ/ssqtot);
model.XRes.EffectLabel='Res';
model.XRes.TermOrder=max(desorder)+1;

% ------------------------------------------------------------
% Type III SSQ and effect matrices (computed on the SELECTED/GROUPED blocks)
% H_obs(block) = (Pfull - Pred_block) * Xp, where Pred_block is the projection
% of the reduced model that excludes the block.
% ------------------------------------------------------------
nBlocks = length(dmatrices);
Dfull  = [];
colIdx = cell(nBlocks,1);
c0 = 0;
for bi = 1:nBlocks
    Di = dmatrices{bi};
    Dfull = [Dfull, Di];
    colIdx{bi} = (c0+1):(c0+size(Di,2));
    c0 = c0 + size(Di,2);
end
Pfull = Dfull * pinv(Dfull);
rank_full = rank(Dfull);

% store minimal bookkeeping (avoid storing huge projection matrices if undesired)
model.TypeIII.ColIdx = colIdx;
model.TypeIII.RankFull = rank_full;

% save block column indices
for bi = 1:nBlocks
    model.Blocks(bi).ColIdx = colIdx{bi};
end

ssq_centered = model.Xdata.CenteredSSQ;

for bi = 2:nBlocks  % usually skip Mean
    keepCols = true(1, size(Dfull,2));
    keepCols(colIdx{bi}) = false;

    Dred = Dfull(:, keepCols);
    Pred = Dred * pinv(Dred);

    Hobs = (Pfull - Pred) * Xp;
    ssqIII = sum(Hobs(:).^2);
    dfIII  = rank_full - rank(Dred);

    l = strtrim(model.TermLabels{bi});

    eval(['model.X', l, '.TypeIII_EffectMatrix = Hobs;']);
    eval(['model.X', l, '.TypeIII_SSQ = ssqIII;']);
    eval(['model.X', l, '.TypeIII_df = dfIII;']);
    eval(['model.X', l, '.TypeIII_ExplVar = 100*(ssqIII/ssq_centered);']);
end



for i=2:length(dmatrices)
    l=strtrim(deslabels{i});

    % Determine which other terms to remove to form the reduced matrix
    if ischar(opt.reducedmodel) && strcmp(opt.reducedmodel, 'standard')
        remLabels = deslabels([2:i-1 i+1:end]);   % all other modeled terms (exclude Mean)
    else
        rem = opt.reducedmodel{i};
        if isempty(rem)
            remLabels = {};
        elseif isnumeric(rem)
            % interpret as indices into current term list
            remLabels = deslabels(rem);
        else
            % interpret as labels (cell/char/string)
            if ischar(rem) || isstring(rem)
                rem = {char(rem)};
            end
            remLabels = cell(size(rem));
            for rr = 1:length(rem)
                s = strtrim(char(rem{rr}));
                % allow grouped specification using '+' or '_' separators
                if contains(s,'+')
                    parts = strsplit(s,'+');
                    parts = cellfun(@canonical_label, parts, 'UniformOutput', false);
                    s = strjoin(parts,'_');
                elseif contains(s,'_')
                    parts = strsplit(s,'_');
                    parts = cellfun(@canonical_label, parts, 'UniformOutput', false);
                    s = strjoin(parts,'_');
                else
                    s = canonical_label(s);
                end
                remLabels{rr} = s;
            end
        end
    end

    Xx = model.Xdata.CenteredData;
    for j=1:length(remLabels)
        m = strtrim(remLabels{j});
        eval(['Xx = Xx - model.X', m, '.EffectMatrix;'])
    end
    eval(['model.X', l, '.ReducedMatrix = Xx;'])
end


model.TermLabels=deslabels;
model.Options=opt;


%Permutation tests, if required
model=ascaptest(model);

%SCA Models
model=ascasca(model);

%bootstrap, if required
model=ascaboot(model);









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dmatrices,desterms, deslabels, desorder]=createdesign(designmat)
nfact=length(designmat);
ns=size(designmat{1},1);


indmat=fullfact(repmat(2,1,nfact));
nmat=size(indmat,1);
dmatrices=cell(1,nmat);
desterms=cell(1,nmat);
deslabels=cell(1,nmat);
desorder=zeros(nmat,1);



for i=1:nmat
    dm=ones(ns,1);
    for j=1:nfact
        if indmat(i,j)==1
            effmat=designmat{j};
        else
            effmat=ones(ns,1);
        end
        dm=kron(dm,effmat);
        dm=dm(1:ns+1:end,:);
    end
    dmatrices{i}=dm;
    desterms{i}=find(indmat(i,:)==1);
    deslabels{i}=char(64+find(indmat(i,:)==1));
    desorder(i)=length(desterms{i});


end

deslabels=sort(char(deslabels),2);


% Canonicalize labels: make BA -> AB, CBA -> ABC, etc.
lab = sort(char(deslabels), 2);   % char matrix
labCell = cellstr(lab);          % cellstr (trim later if you want)

% Sort terms by interaction order (1-way, 2-way, 3-way...) then alphabetically
labChar = char(labCell);         % padded char matrix for consistent width
key = [desorder(:), double(labChar)];
[~, newindex] = sortrows(key);

% Reorder everything consistently
desterms  = desterms(newindex);
dmatrices = dmatrices(newindex);
desorder  = desorder(newindex);
labCell   = labCell(newindex);

% Back to your original label format
deslabels = labCell';
deslabels{1} = 'Mean';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = canonical_label(s)
% Canonicalize a term label (e.g., 'BA' -> 'AB'). Supports single-letter factors.
% Leaves 'Mean' unchanged. For grouped labels like 'C+AC', canonicalizes each member.
if isstring(s); s = char(s); end
s = strtrim(s);
if strcmpi(s,'Mean')
    s = 'Mean';
    return
end
if contains(s,'+')
    parts = strsplit(s,'+');
    parts = cellfun(@canonical_label, parts, 'UniformOutput', false);
    s = strjoin(parts, '+');
    return
end
% keep only letters, sort them
c = char(s);
c = c(c~=' ');
c = sort(c);
s = strtrim(c);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sc=TPstr2cell(s)
n=size(s,1);
sc=cell(n,1);

for i=1:n
    sc{i}=s(i,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dmain=gendmat(designmat)
ns=size(designmat,1);
nfact=size(designmat,2);
dmain=cell(1,nfact);

for i=1:nfact
    lev=unique(designmat(:,i));
    nl=length(lev);
    dmat=zeros(ns,nl-1);
    for j=1:nl-1
        dmat(designmat(:,i)==lev(j),j)=1;
    end
    dmat(designmat(:,i)==lev(nl),:)=-1;
    dmain{i}=dmat;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pmodel=ascaptest(ascamodel)

pmodel = ascamodel;
dlab   = ascamodel.TermLabels;

if strcmp(ascamodel.Options.permtest, 'on')
    signfacts = cell(length(dlab)-1,1);
    sc = 0;
else
    signfacts = [];
end

% --- Build FULL design matrix once (concatenate all term design matrices) ---
Xp = ascamodel.Xdata.PreprData;

Dfull  = [];
colIdx = cell(length(dlab),1);
c0 = 0;

for k = 1:length(dlab)
    lk = strtrim(dlab{k});
    eval(['Dk = ascamodel.X', lk, '.DesignMatrix;']);
    Dfull = [Dfull, Dk];
    colIdx{k} = (c0+1):(c0+size(Dk,2));
    c0 = c0 + size(Dk,2);
end

Pfull = Dfull * pinv(Dfull);

% --- Loop terms (skip Mean) ---
for i = 2:length(dlab)
    l = strtrim(dlab{i});

    if strcmp(ascamodel.Options.permtest, 'on')

        doTest = false;
        if ischar(ascamodel.Options.permfacts) && strcmp(ascamodel.Options.permfacts,'all')
            doTest = true;
        elseif ~ischar(ascamodel.Options.permfacts) && ismember(l, ascamodel.Options.permfacts)
            doTest = true;
        end

        if doTest
            % Reduced design: remove columns of term i, keep all others (incl Mean)
            keep = true(1, size(Dfull,2));
            keep(colIdx{i}) = false;
            Dred = Dfull(:, keep);
            Pred = Dred * pinv(Dred);

            % Observed Type III SSQ for term i:
            Hobs = (Pfull - Pred) * Xp;
            seff = sum(Hobs(:).^2);

            % Permutation null distribution
            if isfield(ascamodel.Options,'permmethod') && strcmpi(ascamodel.Options.permmethod,'freedman-lane')
                ssqp = ptest_fl(Xp, Pfull, Pred, ascamodel.Options.nperm);
            else
                % legacy behavior (not Freedman–Lane)
                ssqp = ptest(Xp, Dfull, colIdx{i}, ascamodel.Options.nperm);
            end

            p = (1 + sum(ssqp >= seff)) ./ (ascamodel.Options.nperm + 1);

            if p <= 0.05
                sc = sc + 1;
                signfacts{sc} = l;
            end

            eval(['pmodel.X', l, '.EffectSignif.NullDistr = ssqp;'])
            eval(['pmodel.X', l, '.EffectSignif.p = p;'])

            % also store Type III SSQ used in the test:
            eval(['pmodel.X', l, '.TypeIII_SSQ = seff;'])

        else
            eval(['pmodel.X', l, '.EffectSignif.NullDistr = [];'])
            eval(['pmodel.X', l, '.EffectSignif.p = [];'])
        end

    else
        eval(['pmodel.X', l, '.EffectSignif.NullDistr = [];'])
        eval(['pmodel.X', l, '.EffectSignif.p = [];'])
    end
end

if strcmp(ascamodel.Options.permtest, 'on')
    signfacts = signfacts(1:sc);
end

pmodel.SignificantTerms = signfacts;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ssqp = ptest_fl(X, Pfull, Pred, nperm)
% Freedman–Lane: permute reduced-model residuals, refit full vs reduced

ns   = size(X,1);
ssqp = zeros(nperm,1);

Xhat_red = Pred * X;
Ered     = X - Xhat_red;

for i = 1:nperm
    hh = randperm(ns);
    Y  = Xhat_red + Ered(hh,:);           % pseudo data under H0
    H  = (Pfull - Pred) * Y;              % Type III effect matrix on Y
    ssqp(i) = sum(H(:).^2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ssqp=ptest(X,D, nperm)

ns=size(X,1);  %Number of samples
ssqp=zeros(nperm,1);   %Initialization of the permuted SSQ vector

for i=1:nperm
    hh=randperm(ns);
    Xpp=D(hh,:)*pinv(D(hh,:))*X;
    ssqp(i)=sum(sum(Xpp.^2));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function smodel=ascasca(ascamodel)

smodel=ascamodel;
dlab=ascamodel.TermLabels;

for i=2:length(dlab)
    l=strtrim(dlab{i});

    if isfield(ascamodel.Options,'sca_source') && strcmpi(ascamodel.Options.sca_source,'typeiii')
        eval(['Xr=ascamodel.X', l, '.TypeIII_EffectMatrix; '])
        eval(['ssqr=ascamodel.X', l, '.TypeIII_SSQ; '])
    else
        eval(['Xr=ascamodel.X', l, '.EffectMatrix; '])
        eval(['ssqr=ascamodel.X', l, '.EffectSSQ; '])
    end

    R=rank(Xr);
    [u,s,P]=svds(Xr,R);
    t=u*s;
    taug=(Xr+ascamodel.XRes.EffectMatrix)*P;
    varex=100*(diag(s).^2)/ssqr;

    eval(['smodel.X', l, '.SCA.Model.Scores=t; '])
    eval(['smodel.X', l, '.SCA.Model.ScoreswithRes=taug; '])
    eval(['smodel.X', l, '.SCA.Model.Loadings=P; '])
    eval(['smodel.X', l, '.SCA.Model.ExplVar=varex; '])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bmodel=ascaboot(ascamodel)

bmodel=ascamodel;
dlab=ascamodel.TermLabels;

useTypeIII = isfield(ascamodel.Options,'sca_source') && strcmpi(ascamodel.Options.sca_source,'typeiii');

if useTypeIII
    % Build FULL design projection once (used for all terms)
    Xp = ascamodel.Xdata.PreprData;

    Dfull  = [];
    colIdx = cell(length(dlab),1);
    c0 = 0;

    for k = 1:length(dlab)
        lk = strtrim(dlab{k});
        eval(['Dk = ascamodel.X', lk, '.DesignMatrix;']);
        Dfull = [Dfull, Dk];
        colIdx{k} = (c0+1):(c0+size(Dk,2));
        c0 = c0 + size(Dk,2);
    end

    Pfull = Dfull * pinv(Dfull);
end

for i=2:length(dlab)
    l=strtrim(dlab{i});

    if strcmp(ascamodel.Options.bootstrap, 'all');
        eval(['Pd=ascamodel.X', l, '.SCA.Model.Loadings; '])

        if useTypeIII
            % Reduced model for this term (remove the term columns)
            keep = true(1, size(Dfull,2));
            keep(colIdx{i}) = false;
            Dred = Dfull(:, keep);
            Pred = Dred * pinv(Dred);

            % Bootstrap around Type III loadings using Freedman–Lane pseudo data
            [Pb, Pbcrit, svars] = bootload_typeiii_fl(Xp, Pfull, Pred, Dred, Pd, ascamodel.Options.confl, ascamodel.Options.nboot);
        else
            % Classic ASCA bootstrap (your existing behavior)
            if strcmp(ascamodel.Options.bootmatrix, 'original')
                Xd = ascamodel.Xdata.CenteredData;
            elseif strcmp(ascamodel.Options.bootmatrix, 'reduced')
                eval(['Xd=ascamodel.X', l, '.ReducedMatrix; '])
            end

            eval(['Dd=ascamodel.X', l, '.DesignMatrix; '])
            [Pb, Pbcrit, svars] = bootload(Xd, Dd, Pd, ascamodel.Options.confl, ascamodel.Options.nboot);
        end

    elseif strcmp(ascamodel.Options.bootstrap, 'signif');
        if ismember(l, ascamodel.SignificantTerms)
            eval(['Pd=ascamodel.X', l, '.SCA.Model.Loadings; '])

            if useTypeIII
                % Reduced model for this term (remove the term columns)
                keep = true(1, size(Dfull,2));
                keep(colIdx{i}) = false;
                Dred = Dfull(:, keep);
                Pred = Dred * pinv(Dred);

                % Bootstrap around Type III loadings using Freedman–Lane pseudo data
                [Pb, Pbcrit, svars] = bootload_typeiii_fl(Xp, Pfull, Pred, Dred, Pd, ascamodel.Options.confl, ascamodel.Options.nboot);
            else
                % Classic ASCA bootstrap (your existing behavior)
                if strcmp(ascamodel.Options.bootmatrix, 'original')
                    Xd = ascamodel.Xdata.CenteredData;
                elseif strcmp(ascamodel.Options.bootmatrix, 'reduced')
                    eval(['Xd=ascamodel.X', l, '.ReducedMatrix; '])
                end

                eval(['Dd=ascamodel.X', l, '.DesignMatrix; '])
                [Pb, Pbcrit, svars] = bootload(Xd, Dd, Pd, ascamodel.Options.confl, ascamodel.Options.nboot);
            end
        else
            Pb=[];
            Pbcrit=[];
            svars=[];
        end


    elseif strcmp(ascamodel.Options.bootstrap, 'off');
        Pb=[];
        Pbcrit=[];
        svars=[];
    end


    switch ascamodel.Options.bootsave
        case 'all'
            eval(['bmodel.X', l, '.SCA.Bootstrap.Loadings=Pb; '])
            eval(['bmodel.X', l, '.SCA.Bootstrap.ConfIntervals=Pbcrit; '])
            eval(['bmodel.X', l, '.SCA.Bootstrap.SignificantVariables=svars; '])
        case 'confint'
            eval(['bmodel.X', l, '.SCA.Bootstrap.Loadings=[]; '])
            eval(['bmodel.X', l, '.SCA.Bootstrap.ConfIntervals=Pbcrit; '])
            eval(['bmodel.X', l, '.SCA.Bootstrap.SignificantVariables=svars; '])
        case 'signvars'
            eval(['bmodel.X', l, '.SCA.Bootstrap.Loadings=[]; '])
            eval(['bmodel.X', l, '.SCA.Bootstrap.ConfIntervals=[]; '])
            eval(['bmodel.X', l, '.SCA.Bootstrap.SignificantVariables=svars; '])
    end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Pb, Pbcrit,svars]=bootload(Xd,Dd, Pd, confl, nboot)
%Bootstrap

sl=(confl+1)/2;
ll=1-sl;
svars=cell(1,size(Pd,2));

Pb=zeros([nboot size(Pd)]);

bootp=zeros(size(Xd,1),1);
lev=unique(Dd, 'rows');

for i=1:nboot
    for j=1:size(lev,1)
        xx=find(ismember(Dd, lev(j,:), 'rows')==1);
        yy=ceil(length(xx)*rand(1,length(xx)));
        bootp(xx)=xx(yy);
    end
    Xboot=Xd(bootp,:);
    Xdp=Dd*pinv(Dd)*Xboot;
    [~,~,vb]=svds(Xdp,size(Pd,2));
    [~,Pb(i,:,:)]=orth_proc(Pd,vb);
end

Pb_sorted = sort(Pb, 1);

lo = max(1, floor(ll*(nboot+1)));
hi = min(nboot, ceil(sl*(nboot+1)));

Pbcrit = squeeze(Pb_sorted([lo hi],:,:));  % 2 x nv x ncomp

for i=1:length(svars)
    svars{i}=find(sign(squeeze(Pbcrit(1,:,i).*Pbcrit(2,:,i)))==1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Pb, Pbcrit, svars] = bootload_typeiii_fl(X, Pfull, Pred, Dred, Pd, confl, nboot)
% Bootstrap loadings around Type III effect matrix using Freedman–Lane style pseudo-data:
%   Y = Xhat_red + Ered*  (Ered* is resampled within reduced-model strata)
%   H* = (Pfull - Pred) * Y
% Then compute loadings of H* and align to Pd with orthogonal Procrustes.

sl = (confl+1)/2;
ll = 1-sl;

ns    = size(X,1);
nv    = size(X,2);
ncomp = size(Pd,2);

Pb = zeros(nboot, nv, ncomp);

% Reduced fit + residuals
Xhat_red = Pred * X;
Ered     = X - Xhat_red;

% Strata defined by reduced design rows (conservative + matches reduced structure)
lev = unique(Dred, 'rows');

for b = 1:nboot
    bootp = zeros(ns,1);

    for j = 1:size(lev,1)
        idx = find(ismember(Dred, lev(j,:), 'rows')==1);
        r   = randi(length(idx), length(idx), 1);   % resample with replacement
        bootp(idx) = idx(r);
    end

    Y = Xhat_red + Ered(bootp,:);     % pseudo data under reduced model
    H = (Pfull - Pred) * Y;           % Type III effect matrix on pseudo data

    % SVD loadings
    k = min(ncomp, rank(H));
    vb = zeros(nv, ncomp);
    if k > 0
        [~,~,vbk] = svds(H, k);
        vb(:,1:k) = vbk;
    end

    % Align to the observed loadings
    [~, vbalign] = orth_proc(Pd, vb);
    Pb(b,:,:) = vbalign;
end

% Empirical CI
Pb_sorted = sort(Pb, 1);
lo = max(1, floor(ll*(nboot+1)));
hi = min(nboot, ceil(sl*(nboot+1)));
Pbcrit = squeeze(Pb_sorted([lo hi],:,:));   % 2 x nv x ncomp

% Significant vars: CI excludes 0
svars = cell(1, ncomp);
for c = 1:ncomp
    svars{c} = find(sign(Pbcrit(1,:,c) .* Pbcrit(2,:,c)) == 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Dcoded, info] = normalize_design_input(Din)
%NORMALIZE_DESIGN_INPUT Normalize design input to consecutive numeric levels.
% Supports:
%   - numeric matrix with non-consecutive level values
%   - cell array (or string array) of labels
% Returns:
%   Dcoded : numeric matrix with levels coded 1..k per factor
%   info   : struct with per-factor mapping and original values

if isstring(Din)
    Din = cellstr(Din);
end

info = struct();
info.InputClass = class(Din);

if isnumeric(Din)
    Dcoded = Din;
    if any(isnan(Din(:)))
        error('asca:DesignMatrixNaN', 'Design matrix contains NaN values. Please remove or impute them before ASCA.');
    end

    [ns, nf] = size(Din);
    info.nSamples = ns;
    info.nFactors = nf;
    info.Factors = repmat(struct('Type','','OriginalLevels',[],'CodeLevels',[]), 1, nf);

    for f = 1:nf
        lev = unique(Din(:,f), 'stable');
        nl  = length(lev);
        Dcoded(:,f) = zeros(ns,1);
        for j = 1:nl
            Dcoded(Din(:,f) == lev(j), f) = j;
        end
        info.Factors(f).Type = 'numeric';
        info.Factors(f).OriginalLevels = lev(:)';
        info.Factors(f).CodeLevels = 1:nl;
    end

elseif iscell(Din)
    % Expect ns-by-nf cell array of labels (strings/chars)
    [ns, nf] = size(Din);
    info.nSamples = ns;
    info.nFactors = nf;
    info.Factors = repmat(struct('Type','','OriginalLevels',{{}},'CodeLevels',[]), 1, nf);

    Dcoded = zeros(ns, nf);

    for f = 1:nf
        col = Din(:,f);

        % Convert any string elements to char for consistent comparison
        for r = 1:ns
            if isstring(col{r})
                col{r} = char(col{r});
            end
        end

        lev = unique(col, 'stable');
        nl  = length(lev);

        for j = 1:nl
            Dcoded(strcmp(col, lev{j}), f) = j;
        end

        info.Factors(f).Type = 'label';
        info.Factors(f).OriginalLevels = lev(:)';
        info.Factors(f).CodeLevels = 1:nl;
    end

else
    error('asca:InvalidDesignInput', 'Design input must be a numeric matrix or a cell array (or string array) of labels.');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r,yrot]=orth_proc(x,y)
%computes orthogonal procrustes rotation projecting matrix y onto the
%subspace spanned by matrix x.
% syntax: [r,yrot]=orth_proc(x,y)
% where r is the rotation matrix and yrot is the procrustes rotated y

[u,~,v]=svd(y'*x, 0);
r=u*v';
yrot=y*r;