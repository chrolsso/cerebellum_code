csv = 'long_df.csv'; % dataframe exported from SimpleAnalysis.m, or python script
measures = {'thickness_median_mm', 'gm_volume_cm3'};%one row of panels per measure
regions = {'all', 'left', 'right'};% one column of panels per region
groups = {'control', 'patient'};
colors = [0.2 0.45 0.8; 0.85 0.3 0.2];%groupcol
label_subjects = true;% write the subject name next to each point
save_as = 'plot.png';% '' = do not save

df = readtable(csv);
df.group = string(df.group); 
df.region = string(df.region);
df.measure = string(df.measure);
df.subject = string(df.subject);

figure('Color', 'w', 'Position', [100 100 330 * numel(regions) 280 * numel(measures)]);


for i = 1:numel(measures)
    for j = 1:numel(regions)
        subplot(numel(measures), numel(regions), (i - 1) * numel(regions) + j);
        hold on
        for g = 1:numel(groups)
            s = df(df.measure == measures{i} & df.region == regions{j} & df.group == groups{g}, :);
            x = g + 0.12 * (rand(height(s), 1) - 0.5);
            scatter(x, s.value, 50, colors(g, :), 'filled');
            if height(s) > 0
                plot(g + [-0.25 0.25], mean(s.value) * [1 1], 'Color', colors(g, :), 'LineWidth', 2);
            end
            if label_subjects
                text(x + 0.05, s.value, extractAfter(s.subject, strlength(s.subject) - 2), 'FontSize', 8);
            end
        end

        set(gca, 'XTick', 1:numel(groups), 'XTickLabel', groups, 'XLim', [0.5 numel(groups) + 0.5], 'Box', 'off');
        title(strrep(regions{j}, '_', ' '));
        if j == 1
            ylabel(strrep(measures{i}, '_', ' '));
        end
    end
end

if ~isempty(save_as)
    exportgraphics(gcf, save_as, 'Resolution', 150); 
end
