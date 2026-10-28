% Plots the files written by 9a and 9b.
outdir = 'D:\SCA_Maurizio_Post\Postprocessed\output';
set(groot, 'defaultFigureColor', 'w');
measure = 'cortex_volume_tiv_pct';       % or thickness_median_mm
label = 'Cortex volume (% of TIV)';
names_file = fullfile(outdir, 'lobule_names.csv');     % optional, columns: label,name

R = readtable(fullfile(outdir, 'stats', 'nighres_patient_vs_control.csv'));
R.region = string(R.region); R.measure = string(R.measure);

% 1. difference per lobule, patient minus control, adjusted for age (% of control mean, 95% CI)
F = R(R.measure == measure & startsWith(R.region, 'lobule_'), :);
F.rel_low = 100 * F.ci_low ./ F.control_mean; F.rel_high = 100 * F.ci_high ./ F.control_mean;
F = sortrows(F, 'rel_diff_pct');
lab = F.region;
if isfile(names_file)
    N = readtable(names_file); num = str2double(extractAfter(F.region, 'lobule_'));
    for k = 1:height(F), j = find(N.label == num(k), 1); if ~isempty(j), lab(k) = string(N.name(j)); end, end
end
figure(901); clf; hold on
for k = 1:height(F)
    plot([F.rel_low(k) F.rel_high(k)], [k k], 'k', 'LineWidth', 1.5);
    if F.p_fdr(k) < 0.05, fc = 'k'; else, fc = 'w'; end
    plot(F.rel_diff_pct(k), k, 'ko', 'MarkerFaceColor', fc);
end
xline(0, '--'); set(gca, 'YTick', 1:height(F), 'YTickLabel', lab, 'YDir', 'reverse'); box off
xlabel('Patient minus control (% of control mean), filled = FDR < 0.05'); title(label)

% 2. whole, left and right cerebellum, adjusted for age (same as 8c)
Lg = readtable(fullfile(outdir, 'nighres_long.csv'));
Lg.region = string(Lg.region); Lg.measure = string(Lg.measure); Lg.status = string(Lg.status);
raw = strrep(measure, '_tiv_pct', '_cm3'); regs = {'all', 'left', 'right'};
figure(902); clf; col = lines(2);
for r = 1:numel(regs)
    S = Lg(Lg.measure == raw & Lg.region == regs{r}, :);
    y = S.value; if endsWith(measure, '_tiv_pct'), y = 100 * y ./ S.tiv; end
    D = rmmissing(table(y, S.age, double(S.status == "Patient"), 'VariableNames', {'y', 'age', 'patient'}));
    m = fitlm(D, 'y ~ age + patient');
    D.adj = D.y - m.Coefficients{'age', 'Estimate'} * (D.age - mean(D.age));
    subplot(1, numel(regs), r); hold on
    for g = 1:2
        yy = D.adj(D.patient == g - 1); mu = mean(yy); sd = std(yy);
        scatter(g + 0.15 * (rand(size(yy)) - 0.5) * 2, yy, 45, col(g, :), 'filled', 'MarkerFaceAlpha', 0.7);
        plot(g + [-0.25 0.25], [mu mu], 'k', 'LineWidth', 2); plot([g g], [mu - sd mu + sd], 'k', 'LineWidth', 1.5);
    end
    set(gca, 'XTick', 1:2, 'XTickLabel', {'Control', 'Patient'}, 'XLim', [0.4 2.6]); box off
    title(sprintf('%s (p=%.3f)', regs{r}, m.Coefficients{'patient', 'pValue'}))
    if r == 1, ylabel([label ', adjusted to mean age']); end
end

% 3. nighres against FastSurfer cortex volume (check, not a result)
B = readtable(fullfile(outdir, 'stats', 'nighres_vs_fastsurfer.csv'));
figure(903); clf
gscatter(B.fastsurfer_cortex_cm3, B.nighres_cortex_cm3, B.status)
xlabel('FastSurfer cerebellar cortex (cm^3)'); ylabel('Nighres cerebellar cortex (cm^3)'); box off
