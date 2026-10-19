function [hPts, hEll, out] = plotasca_bootstrapellipses(factormodel, sccomp, confLevel, B, whichScores)
% plotasca_bootstrapellipses Plot ASCA scores with bootstrap mean ellipses
%
% Bootstrap procedure (per group):
%   for b=1..B: resample rows within group, compute mean score mu*(b)
% Ellipse: (z-muhat)' inv(Cov_boot(mu*)) (z-muhat) <= chi2inv(confLevel,2)

if nargin < 2 || isempty(sccomp), sccomp = [1 2]; end
if nargin < 3 || isempty(confLevel), confLevel = 0.95; end
if nargin < 4 || isempty(B), B = 2000; end
if nargin < 5 || isempty(whichScores), whichScores = 'Scores'; end

T = factormodel.SCA.Model.(whichScores);
T2 = T(:, sccomp(1:2));

[~,~,lindex] = unique(factormodel.DesignMatrix,'rows','stable');
ng = max(lindex);

Mark = {'o','s','d','^','v','>','<','p','h','+','*','.','x'};
colr = lines(7);

figure; hold on;
hPts = gobjects(ng,1);
hEll = gobjects(ng,1);

theta = linspace(0,2*pi,200);
unitCircle = [cos(theta); sin(theta)];
c = chi2inv(confLevel, 2);

out = struct();
out.confLevel = confLevel;
out.B = B;
out.groups = cell(ng,1);

for g = 1:ng
    idx = (lindex==g);
    Y = T2(idx,:); n = size(Y,1);

    f1 = rem(g,13); if f1==0, f1=13; end
    f2 = rem(g,7);  if f2==0, f2=7;  end

    hPts(g) = plot(Y(:,1), Y(:,2), Mark{f1}, ...
        'markerfacecolor', colr(f2,:), 'markeredgecolor', colr(f2,:), 'markersize', 10);

    if n < 2
        out.groups{g} = struct('n',n,'mu',mean(Y,1),'CovBoot',nan(2),'c',c,'ellipse',nan(2,numel(theta)));
        continue;
    end

    muhat = mean(Y,1);

    % Bootstrap distribution of group mean
    muBoot = zeros(B,2);
    for b = 1:B
        ii = randi(n, [n,1]);          % resample indices with replacement
        muBoot(b,:) = mean(Y(ii,:),1);
    end

    CovBoot = cov(muBoot,0);          % covariance of bootstrap means
    CovBoot = (CovBoot + CovBoot')/2;

    % Ellipse from CovBoot
    [R,flag] = chol(CovBoot,'lower');
    if flag ~= 0
        [V,D] = eig(CovBoot);
        D = max(D,0);
        R = V*sqrt(D);
    end

    ell = muhat(:) + sqrt(c) * (R * unitCircle);
    hEll(g) = plot(ell(1,:), ell(2,:), '-', 'linewidth', 1.8);

    out.groups{g} = struct('n',n,'mu',muhat,'CovBoot',CovBoot,'c',c,'ellipse',ell,'muBoot',muBoot);
end

xx = axis;
plot(xx(1:2), [0 0], '--k', 'linewidth', 1.2);
plot([0 0], xx(3:4), '--k', 'linewidth', 1.2);

xlabel(['Scores on SC', num2str(sccomp(1))], 'fontsize', 12, 'fontweight', 'bold');
ylabel(['Scores on SC', num2str(sccomp(2))], 'fontsize', 12, 'fontweight', 'bold');
set(gca, 'linewidth', 1.5, 'fontsize', 12, 'fontweight', 'bold');

title(sprintf('ASCA %s with bootstrap %.1f%% mean ellipses (B=%d)', whichScores, 100*confLevel, B), ...
      'fontsize', 12, 'fontweight', 'bold');

end
