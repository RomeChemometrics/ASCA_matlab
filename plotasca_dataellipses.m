function [hPts, hEll, out] = plotasca_dataellipses(factormodel, sccomp, coverage, whichScores)
% plotasca_dataellipses Plot ASCA scores with data ellipses per group
%
% Data ellipse: (z-mu)' * inv(S) * (z-mu) <= chi2inv(coverage,2)

if nargin < 2 || isempty(sccomp), sccomp = [1 2]; end
if nargin < 3 || isempty(coverage), coverage = 0.95; end
if nargin < 4 || isempty(whichScores), whichScores = 'Scores'; end

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
c = chi2inv(coverage, 2);

out = struct();
out.coverage = coverage;
out.groups = cell(ng,1);

for g = 1:ng
    idx = (lindex==g);
    Y = T2(idx,:); n = size(Y,1);

    f1 = rem(g,13); if f1==0, f1=13; end
    f2 = rem(g,7);  if f2==0, f2=7;  end

    hPts(g) = plot(Y(:,1), Y(:,2), Mark{f1}, ...
        'markerfacecolor', colr(f2,:), 'markeredgecolor', colr(f2,:), 'markersize', 10);

    if n < 3
        out.groups{g} = struct('n',n,'mu',mean(Y,1),'S',nan(2),'c',c,'ellipse',nan(2,numel(theta)));
        continue;
    end

    mu = mean(Y,1);
    S  = cov(Y,0);               % unbiased covariance
    S  = (S+S')/2;

    [R,flag] = chol(S,'lower');
    if flag ~= 0
        [V,D] = eig(S);
        D = max(D,0);
        R = V*sqrt(D);
    end

    ell = mu(:) + sqrt(c) * (R * unitCircle);
    hEll(g) = plot(ell(1,:), ell(2,:), '-', 'linewidth', 1.8);

    out.groups{g} = struct('n',n,'mu',mu,'S',S,'c',c,'ellipse',ell);
end

xx = axis;
plot(xx(1:2), [0 0], '--k', 'linewidth', 1.2);
plot([0 0], xx(3:4), '--k', 'linewidth', 1.2);

xlabel(['Scores on SC', num2str(sccomp(1))], 'fontsize', 12, 'fontweight', 'bold');
ylabel(['Scores on SC', num2str(sccomp(2))], 'fontsize', 12, 'fontweight', 'bold');
set(gca, 'linewidth', 1.5, 'fontsize', 12, 'fontweight', 'bold');

title(sprintf('ASCA %s with %.1f%% data ellipses', whichScores, 100*coverage), ...
      'fontsize', 12, 'fontweight', 'bold');

end
