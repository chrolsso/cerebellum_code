% Plots the files written by 8a and 8b. Core MATLAB only.
merged = 'clinical_merged.csv';          % from 8a
statsdir = 'stats/';                     % from 8b
outcomes = {'cerebellum_total', 'brain_total'};   % one panel each
model = 'cerebellum_by_group_brain';     % model shown in the coefficient plot

T = readtable(merged); T.group = string(T.group);
groups = unique(T.group); if any(groups == "Control"), groups = [groups(groups == "Control"); groups(groups ~= "Control")]; end
col = lines(numel(groups));

% 1. volume per group (box = median and quartiles, dots = subjects)
f1 = figure('Position', [60 60 420 * numel(outcomes) 420]);
for k = 1:numel(outcomes)
    subplot(1, numel(outcomes), k); hold on
    for g = 1:numel(groups)
        y = T.(outcomes{k})(T.group == groups(g)); y = y(~isnan(y));
        x = g + 0.12 * (rand(size(y)) - 0.5) * 2;
        scatter(x, y, 40, col(g, :), 'filled', 'MarkerFaceAlpha', 0.7);
        if numel(y) > 1, q = prctile(y, [25 50 75]); plot(g + [-0.25 0.25], [q(2) q(2)], 'k', 'LineWidth', 2); plot([g g], [q(1) q(3)], 'k'); end
    end
    set(gca, 'XTick', 1:numel(groups), 'XTickLabel', groups, 'XTickLabelRotation', 30); ylabel([strrep(outcomes{k}, '_', ' ') ' (cm^3)']); box off
end
saveas(f1, 'plot_volume_by_group.png');

% 2. cerebellum against brain volume, with a line per group of 3 or more
f2 = figure('Position', [60 60 520 440]); hold on
for g = 1:numel(groups)
    i = T.group == groups(g); x = T.brain_total(i); y = T.cerebellum_total(i);
    scatter(x, y, 50, col(g, :), 'filled', 'DisplayName', char(groups(g)));
    if sum(i) >= 3, p = polyfit(x, y, 1); xs = linspace(min(x), max(x), 20); plot(xs, polyval(p, xs), 'Color', col(g, :), 'HandleVisibility', 'off'); end
end
xlabel('brain volume (cm^3)'); ylabel('cerebellar volume (cm^3)'); legend('Location', 'best'); box off
saveas(f2, 'plot_cerebellum_vs_brain.png');

% 3. model coefficients with 95% CI (group versus control)
C = readtable([statsdir 'coefficients.csv']); C = C(strcmp(C.model, model) & startsWith(C.term, 'C(group'), :);
if ~isempty(C)
    f3 = figure('Position', [60 60 560 120 + 40 * height(C)]); hold on
    for k = 1:height(C)
        plot([C.ci_low(k) C.ci_high(k)], [k k], 'k', 'LineWidth', 1.5); plot(C.estimate(k), k, 'ko', 'MarkerFaceColor', 'k');
    end
    xline(0, '--'); set(gca, 'YTick', 1:height(C), 'YTickLabel', regexprep(C.term, '.*\[T\.(.*)\]', '$1'), 'YDir', 'reverse');
    xlabel('difference from control (cm^3), 95% CI'); title(strrep(model, '_', ' ')); box off
    saveas(f3, 'plot_group_effects.png');
end
