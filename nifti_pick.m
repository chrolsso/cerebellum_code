function files = nifti_pick(startdir)
% NIFTI_PICK  File dialog for choosing NIfTI images, also from several folders: after each selection you can add more
% from another folder. Returns a cell array of full paths ({} if cancelled).
%   files = nifti_pick;   files = nifti_pick('C:\data');   nifti_view(nifti_pick)
if nargin < 1 || isempty(startdir) || ~exist(startdir, 'dir'), startdir = pwd; end
files = {};
while true
    [f, p] = uigetfile({'*.nii;*.nii.gz;*.mgz;*.mgh', 'NIfTI and FreeSurfer images (*.nii, *.nii.gz, *.mgz, *.mgh)'; '*.*', 'All files'}, ...
        sprintf('Select image(s)  [%d chosen so far]', numel(files)), startdir, 'MultiSelect', 'on');
    if isequal(f, 0), break; end
    if ischar(f), f = {f}; end
    files = [files, fullfile(p, f)]; %#ok<AGROW>
    startdir = p;
    answer = questdlg(sprintf('%d image(s) chosen.', numel(files)), 'nifti_pick', 'Add from another folder', 'Done', 'Done');
    if ~strcmp(answer, 'Add from another folder'), break; end
end
end
