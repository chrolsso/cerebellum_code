root = '/';  %root directory, typically Postprocessed/Input
subjects = {'SUBID'};                                          % subject folders
groups = {'patient', 'control'}; %Don't know this one yet, but it should match subjects
%TODO: automate so that it catchers


rows = {};
for s = 1:numel(subjects)
    subject = subjects{s}; 
    group = groups{s};
    [thick, zoom] = read_nifti(find_one(root, subject, 'surfaces', '*cruise-thick.nii.gz'));
    cortex = read_nifti(find_one(root, subject, 'surfaces', '*cruise-cortex.nii.gz')) > 0.5;
    gm = read_nifti(find_one(root, subject, 'tissues', '*_msk_gmrawprior.nii.gz'));
    wm = read_nifti(find_one(root, subject, 'tissues', '*_msk_wmrawprior.nii.gz'));
    voxel_mm = mean(zoom); 
    voxel_cm3 = prod(zoom) / 1000;
    scale = voxel_mm; %OBS check that this is correct!

    regions = {'all',   read_nifti(find_one(root, subject, 'coreg', '*_cr_cbmask_def-img.nii.gz')) > 0.5; ...
               'left',  read_nifti(find_one(root, subject, 'coreg', '*_cr_cbmask_l_def-img.nii.gz')) > 0.5; ...
               'right', read_nifti(find_one(root, subject, 'coreg', '*_cr_cbmask_r_def-img.nii.gz')) > 0.5};

    for r = 1:size(regions, 1)
        region = regions{r, 1}; mask = regions{r, 2};
        t = double(thick(cortex & mask & (thick > 0))) * scale;
        values = {'thickness_median_mm', median(t); 'thickness_mean_mm', mean(t); 'thickness_p10_mm', prctile(t, 10); 'thickness_p90_mm', prctile(t, 90); ...
                  'cortex_volume_cm3', sum(cortex(:) & mask(:)) * voxel_cm3; 'gm_volume_cm3', sum(double(gm(:)) .* mask(:)) * voxel_cm3; ...
                  'wm_volume_cm3', sum(double(wm(:)) .* mask(:)) * voxel_cm3; 'region_volume_cm3', sum(mask(:)) * voxel_cm3};
        for v = 1:size(values, 1)
            rows(end + 1, :) = {subject, group, region, values{v, 1}, values{v, 2}};
        end
    end
end

long_df = cell2table(rows, 'VariableNames', {'subject', 'group', 'region', 'measure', 'value'});
writetable(long_df, [root 'long_df.csv']);
disp(unstack(long_df(long_df.region == "all" | strcmp(long_df.region, 'all'), {'measure', 'group', 'value'}), 'value', 'group', 'AggregationFunction', @mean))


function file = find_one(root, subject, folder, pattern)
f = dir([root subject '/nighres/' folder '/' pattern]);
assert(numel(f) == 1, 'Expected one file: %s %s (found %d)', subject, pattern, numel(f));
file = fullfile(f.folder, f.name);
end

function [data, zoom] = read_nifti(file) 
data=niftiread(file);
info=inftiinfo(file);
zoom=info.PixelDimenstions(1:3);
end
