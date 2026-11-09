% Plots the files written by 10a (profiles per region and a sample of single-voxel profiles).
outdir = 'D:\SCA_Maurizio_Post\Postprocessed\output\laminar';
subject_id = 'ID';
set(groot, 'defaultFigureColor', 'w');
region = 'all';                  % 'all' or e.g. 'lobule_12'
n_show = 300;                    % single-voxel profiles to draw
layer_to_show = 4;               % sampled surface to map in the slice view (0 = WM side, last = CSF side)
compare = {'lobule_3', 'lobule_5', 'lobule_12'};

P = readtable(fullfile(outdir, [subject_id '_laminar_profile.csv']));      % mean, sd and gradient per region and depth
V = readtable(fullfile(outdir, [subject_id '_laminar_voxels.csv']));       % sample of single-voxel profiles
P.region = string(P.region); V.region = string(V.region);
cols = V.Properties.VariableNames(startsWith(V.Properties.VariableNames, 't1_'));
depth = linspace(0, 1, numel(cols));
getv = @(r) V{V.region == string(r), cols};                                % voxels x layers
X = getv(region); fprintf('%s: %d voxel profiles\n', region, size(X, 1));

% 1. individual voxel profiles with mean and median
figure(1001); clf; hold on
pick = X(randperm(size(X, 1), min(n_show, size(X, 1))), :);
plot(depth, pick', 'Color', [0.6 0.6 0.6 0.2]);
plot(depth, mean(X), 'r', 'LineWidth', 2); plot(depth, median(X), 'b--', 'LineWidth', 2);
xlabel('relative depth (0 = WM side, 1 = CSF side)'); ylabel('T1'); title(['Voxel profiles, ' char(region)]); legend({'voxels', 'mean', 'median'}); box off

% 2. distribution of the intensities at each depth (density with 10/50/90 % lines)
figure(1002); clf
edges = linspace(prctile(X(:), 0.5), prctile(X(:), 99.5), 60);
D = zeros(numel(edges) - 1, numel(depth));
for i = 1:numel(depth), D(:, i) = histcounts(X(:, i), edges, 'Normalization', 'pdf'); end
imagesc(depth, edges(1:end-1) + diff(edges) / 2, D); axis xy; colormap(parula); colorbar; hold on
plot(depth, prctile(X, 50), 'w', 'LineWidth', 1.5); plot(depth, prctile(X, 10), 'w:', 'LineWidth', 1.2); plot(depth, prctile(X, 90), 'w:', 'LineWidth', 1.2);
xlabel('relative depth'); ylabel('T1'); title(['T1 distribution per depth (white: 10/50/90 %), ' char(region)])

% 3. box plot at every sampled depth
figure(1003); clf
boxplot(X, 'Positions', depth, 'Widths', 0.04, 'Symbol', ''); set(gca, 'XTick', depth, 'XTickLabel', compose('%.2f', depth)); box off
xlabel('relative depth'); ylabel('T1'); title(['T1 per sampled surface, ' char(region)])

% 4. where on the cortex: voxels of the sample in one coronal slice, coloured by the T1 on the chosen surface
figure(1004); clf
Vr = V(V.region == string(region), :);
slice_j = round(median(Vr.j)); near = abs(Vr.j - slice_j) <= 3;        % a thin slab of the sample around the median slice
scatter(Vr.i(near), Vr.k(near), 8, Vr.(cols{layer_to_show + 1})(near), 'filled'); axis equal tight; colormap(parula); colorbar
set(gca, 'XTick', [], 'YTick', []); title(sprintf('T1 on surface %d, voxels within 3 slices of j = %d, %s', layer_to_show, slice_j, region))

% 5. compare regions: median profile with a 10-90 % band
figure(1005); clf; hold on; col = lines(numel(compare));
for r = 1:numel(compare)
    Xr = getv(compare{r}); if size(Xr, 1) < 50, fprintf('too few voxels in %s\n', compare{r}); continue; end
    lo = prctile(Xr, 10); hi = prctile(Xr, 90);
    fill([depth fliplr(depth)], [lo fliplr(hi)], col(r, :), 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    plot(depth, median(Xr), 'Color', col(r, :), 'LineWidth', 1.8);
end
xlabel('relative depth'); ylabel('T1'); title('Median profile and 10-90 % band'); box off

% 6. mean profile and its gradient for the whole cortex and each lobule (from the summary file)
figure(1006); clf
regs = unique(P.region, 'stable');
subplot(1, 2, 1); hold on
for r = 1:numel(regs), S = P(P.region == regs(r), :); plot(S.depth, S.mean_t1 / mean(S.mean_t1), 'LineWidth', 0.8 + 1.5 * (regs(r) == "all")); end
xlabel('relative depth'); title('Mean profile (each divided by its own mean)'); box off
subplot(1, 2, 2); hold on
for r = 1:numel(regs), S = P(P.region == regs(r), :); plot(S.depth, S.gradient, 'LineWidth', 0.8 + 1.5 * (regs(r) == "all")); end
yline(0, '--'); xlabel('relative depth'); title('Gradient of the mean profile'); box off
