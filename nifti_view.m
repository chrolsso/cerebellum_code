function fig = nifti_view(files, varargin)
% NIFTI_VIEW  Orthogonal NIfTI viewer modelled on FSLeyes. Core MATLAB only: no SPM, no toolboxes, no downloads.
%
%   nifti_view                          opens a file dialog (nifti_pick): choose images, also from several folders
%   nifti_view('T1.nii.gz')
%   nifti_view({'T1.nii.gz','csfprior.nii.gz'}, 'Colormaps',{'gray','red-yellow'}, 'Ranges',[0.8 5; 0 1])
%   nifti_view(files, 'Snapshot','qc.png', 'Cursor',[260 180 150])      % no window, writes a PNG
%
% Layers: the first file is the base, the others are drawn on top (the list order is the drawing order).
% In the window: "+ add" opens the file dialog again (several folders possible), "- del" removes the selected layer,
% "^" / "v" move it up / down. Every layer has its own colormap, display range, opacity, show and clip-below-min setting.
%
% Name-value options
%   Colormaps    cell, one per file: gray hot jet cool red-yellow blue-lightblue green red blue (default gray, then hot)
%   Ranges       n-by-2 display ranges (default: base 2-99.8 percentile of non-zero voxels; overlays 0-1 if they look like probabilities)
%   Opacity      vector, one per file (default 1)
%   Cursor       voxel [i j k] (1-based) in the FIRST file's own voxel order
%   World        cursor as [x y z] mm instead
%   Space        'auto' (same matrix size -> voxel space, otherwise world), 'voxel' or 'world'.
%                Some pipeline outputs have a slab header although they are on a cropped grid; overlays on the same grid
%                line up in voxel space, which is what 'auto' uses.
%   Volume       which volume of a 4D file to start with (default 1); in the window "," and "." step the selected layer
%   Radiological left-right flipped as in FSLeyes' default (default true)
%   Snapshot     file name: build the figure invisibly, save a PNG and return
%   Hidden       true: build the figure invisibly and return its handle (for scripts and tests)
%
% Mouse: click = move cursor. Wheel / arrow up-down = slice of the active view (left-right: 10 slices; page up/down too).
% Keys: , and . = previous / next volume of the selected layer (4D files), [ and ] = opacity of the selected layer -/+ 0.1, g = radiological/neurological, r = auto range, s = save PNG.

if nargin < 1 || isempty(files), files = nifti_pick; if isempty(files), return; end, end
if ischar(files) || isstring(files), files = cellstr(files); end
ip = inputParser;
addParameter(ip, 'Colormaps', {}); addParameter(ip, 'Ranges', []); addParameter(ip, 'Opacity', []);
addParameter(ip, 'Cursor', []); addParameter(ip, 'World', []); addParameter(ip, 'Snapshot', '');
addParameter(ip, 'Space', 'auto'); addParameter(ip, 'Volume', 1); addParameter(ip, 'Radiological', true);
addParameter(ip, 'Hidden', false);
parse(ip, varargin{:}); o = ip.Results;

tmp = tempname; cleaner = onCleanup(@() cleanup_tmp(tmp)); %#ok<NASGU>

% base image defines the grid and the orientation
[base, Ab, pix, nv1] = read_nifti(files{1}, o.Volume, tmp);
[base, Ab] = deoblique(base, Ab);              % tilted slabs are shown in world space, as nilearn and FSLeyes do
lfile = files; lnv = ones(1, numel(files)); lnv(1) = nv1; lvol = ones(1, numel(files));   % per layer: file, number of volumes, shown volume
sz = size(base);
M = Ab(1:3, 1:3); [~, wax] = max(abs(M), [], 1);
if numel(unique(wax)) < 3, wax = 1:3; end
sgn = ones(1, 3); for c = 1:3, sgn(c) = sign(M(wax(c), c)); if sgn(c) == 0, sgn(c) = 1; end, end
perm = zeros(1, 3); for c = 1:3, perm(wax(c)) = c; end
P = zeros(4); P(4, 4) = 1;
for oo = 1:3
    c = perm(oo);
    if sgn(c) > 0, P(c, oo) = 1; else, P(c, oo) = -1; P(c, 4) = sz(c) - 1; end
end
Aras = Ab * P; vox = abs([Aras(1, 1) Aras(2, 2) Aras(3, 3)]); if any(vox == 0), vox = pix(:)'; end

lvol(1) = min(o.Volume, lnv(1)); vols = {orient(base)}; names = {lname(name_of(files{1}), lvol(1), lnv(1))};
for q = 2:numel(files)
    [v, nm, lnv(q)] = place(files{q}, o.Volume); lvol(q) = min(o.Volume, lnv(q)); vols{q} = v; names{q} = lname(nm, lvol(q), lnv(q)); %#ok<AGROW>
end
n = numel(vols); dim = size(vols{1}); if numel(dim) < 3, dim(3) = 1; end

% layer settings
cycle = {'red-yellow', 'blue-lightblue', 'green', 'hot', 'cool'};
newlayer = @(nm, cm, lo, hi, al, cl) struct('name', nm, 'cmap', cm, 'lo', lo, 'hi', hi, 'alpha', al, 'vis', true, 'clip', cl);
lyr = repmat(newlayer('', 'gray', 0, 1, 1, false), 1, n);
for q = 1:n
    cm = 'gray';
    if numel(o.Colormaps) >= q && ~isempty(o.Colormaps{q}), cm = o.Colormaps{q}; end
    if size(o.Ranges, 1) >= q && ~any(isnan(o.Ranges(q, :))), lo = o.Ranges(q, 1); hi = o.Ranges(q, 2);
    else, [lo, hi] = auto_range(vols{q}, q == 1); end
    al = 1; if numel(o.Opacity) >= q, al = o.Opacity(q); end
    lyr(q) = newlayer(names{q}, cm, lo, hi, al, q > 1);
end

% cursor
if ~isempty(o.World)
    f = Aras \ [o.World(:); 1]; cur = round(f(1:3)' + 1);
elseif ~isempty(o.Cursor)
    f = P \ [o.Cursor(:) - 1; 1]; cur = round(f(1:3)' + 1);
else
    nz = find(vols{1} ~= 0 & ~isnan(vols{1}));
    if isempty(nz), cur = round(dim / 2); else
        [a, b, c] = ind2sub(dim, nz(1:max(1, floor(numel(nz) / 20000)):end)); cur = round([(min(a) + max(a)) / 2, (min(b) + max(b)) / 2, (min(c) + max(c)) / 2]);
    end
end
cur = max(1, min(dim, cur)); sel = 1; act = 3; radio = o.Radiological; lastdir = fileparts(files{1});

% figure
fig = figure('Name', 'nifti_view', 'NumberTitle', 'off', 'Color', 'k', 'MenuBar', 'none', 'ToolBar', 'none', ...
    'Position', [60 60 1280 780]);
if ~isempty(o.Snapshot) || o.Hidden, set(fig, 'Visible', 'off'); end
pos = {[0.02 0.34 0.31 0.62], [0.345 0.34 0.31 0.62], [0.67 0.34 0.31 0.62]};
ax = gobjects(1, 3); him = gobjects(1, 3); hl = gobjects(2, 3); hlab = gobjects(3, 4);
for v = 1:3
    ax(v) = axes('Parent', fig, 'Position', pos{v}, 'Color', 'k', 'XColor', 'none', 'YColor', 'none');
    hold(ax(v), 'on');
    him(v) = image(ax(v), 'CData', zeros(2, 2, 3));
    set(him(v), 'ButtonDownFcn', @(~, ~) click(v)); set(ax(v), 'ButtonDownFcn', @(~, ~) click(v));
    hl(1, v) = line(ax(v), [0 1], [0 1], 'Color', [0 1 0], 'HitTest', 'off');
    hl(2, v) = line(ax(v), [0 1], [0 1], 'Color', [0 1 0], 'HitTest', 'off');
    set(ax(v), 'YDir', 'normal'); daspect(ax(v), [1 1 1]);
    hlab(v, 1) = text(ax(v), 0.01, 0.5, '', 'Units', 'normalized', 'Color', [1 1 .3], 'FontSize', 13, 'HitTest', 'off');
    hlab(v, 2) = text(ax(v), 0.97, 0.5, '', 'Units', 'normalized', 'Color', [1 1 .3], 'FontSize', 13, 'HitTest', 'off');
    hlab(v, 3) = text(ax(v), 0.5, 0.97, '', 'Units', 'normalized', 'Color', [1 1 .3], 'FontSize', 13, 'HitTest', 'off');
    hlab(v, 4) = text(ax(v), 0.5, 0.03, '', 'Units', 'normalized', 'Color', [1 1 .3], 'FontSize', 13, 'HitTest', 'off');
end
bg = [0.15 0.15 0.15]; fg = [0.95 0.95 0.95];
mk = @(style, p, str, cb) uicontrol(fig, 'Style', style, 'Units', 'normalized', 'Position', p, 'String', str, ...
    'Callback', cb, 'BackgroundColor', bg, 'ForegroundColor', fg, 'FontSize', 10);
lab = @(p, str) uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', p, 'String', str, ...
    'BackgroundColor', 'k', 'ForegroundColor', [0.8 0.8 0.8], 'HorizontalAlignment', 'left', 'FontSize', 10);
hlist = mk('listbox', [0.02 0.085 0.23 0.215], list_strings(), @(~, ~) pick());
mk('pushbutton', [0.02 0.02 0.07 0.05], '+ add', @(~, ~) addcb()); mk('pushbutton', [0.095 0.02 0.07 0.05], '- del', @(~, ~) delcb());
mk('pushbutton', [0.17 0.02 0.035 0.05], '^', @(~, ~) movecb(-1)); mk('pushbutton', [0.21 0.02 0.04 0.05], 'v', @(~, ~) movecb(1));
hcmap = mk('popupmenu', [0.27 0.255 0.12 0.04], {'gray', 'hot', 'jet', 'cool', 'red-yellow', 'blue-lightblue', 'green', 'red', 'blue'}, @(~, ~) cmapcb());
hrad = mk('checkbox', [0.41 0.255 0.09 0.04], 'radiol.', @(~, ~) radcb()); set(hrad, 'Value', radio);
lab([0.27 0.195 0.06 0.035], 'min / max'); hlo = mk('edit', [0.335 0.19 0.055 0.045], '', @(~, ~) rangecb()); hhi = mk('edit', [0.395 0.19 0.055 0.045], '', @(~, ~) rangecb());
lab([0.27 0.135 0.06 0.035], 'opacity'); halpha = mk('slider', [0.335 0.135 0.11 0.04], '', @(~, ~) alphacb()); set(halpha, 'Min', 0, 'Max', 1, 'SliderStep', [0.05 0.2]);
halphaval = mk('edit', [0.45 0.13 0.05 0.045], '', @(~, ~) alphaeditcb());
hvis = mk('checkbox', [0.27 0.075 0.06 0.04], 'show', @(~, ~) viscb()); hclip = mk('checkbox', [0.335 0.075 0.11 0.04], 'hide <= min', @(~, ~) clipcb());
mk('pushbutton', [0.27 0.015 0.06 0.045], 'auto', @(~, ~) autocb()); mk('pushbutton', [0.335 0.015 0.06 0.045], 'PNG', @(~, ~) savepng());
hinfo = uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', [0.52 0.02 0.46 0.28], 'String', '', ...
    'BackgroundColor', 'k', 'ForegroundColor', [0.8 1 0.8], 'HorizontalAlignment', 'left', 'FontName', 'Courier New', 'FontSize', 10);
set(fig, 'KeyPressFcn', @key, 'WindowScrollWheelFcn', @wheel, 'WindowButtonMotionFcn', @hover);
setappdata(fig, 'addLayer', @add_layer); setappdata(fig, 'selectLayer', @(k) select_layer(k));
refresh_all();

if ~isempty(o.Snapshot)
    set(fig, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off'); print(fig, o.Snapshot, '-dpng', '-r100'); close(fig); fig = []; return
end

    function v = orient(v)
        v = permute(v, perm);
        for k = 1:3, if sgn(perm(k)) < 0, v = flip(v, k); end, end
    end
    function nm = name_of(file)
        [~, nm, ext] = fileparts(file); if strcmp(ext, '.gz'), [~, nm] = fileparts(nm); end
    end
    function nm = lname(nm, vol, nv)
        if nv > 1, nm = sprintf('%s [%d/%d]', nm, vol, nv); end
    end
    function [v, nm, nvq] = place(file, vol)
        [v, a, ~, nvq] = read_nifti(file, vol, tmp); nm = name_of(file);
        sameSize = isequal(size(v), sz); sameAff = max(abs(a(:) - Ab(:))) < 1e-3;
        mode = o.Space; if strcmp(mode, 'auto'), if sameSize, mode = 'voxel'; else, mode = 'world'; end, end
        if strcmp(mode, 'voxel')
            if ~sameSize, error('Space voxel: %s has a different matrix size than the base image.', nm); end
            if ~sameAff, fprintf('note: %s has another header than the base image; shown in voxel space.\n', nm); end
        elseif ~(sameSize && sameAff)
            fprintf('resampling %s onto the base image grid (world space)...\n', nm);
            v = resample_to(v, a, sz, Ab);
        end
        v = orient(v);
    end
    function add_layer(file)
        [v, nm, nvq] = place(file, o.Volume); n = n + 1; vols{n} = v;
        lfile{n} = file; lnv(n) = nvq; lvol(n) = min(o.Volume, nvq);
        [lo, hi] = auto_range(v, false);
        lyr(n) = newlayer(lname(nm, lvol(n), nvq), 'gray', lo, hi, 1, true); sel = n; lastdir = fileparts(file); refresh_all();
    end
    function setvol(d)                                  % previous / next volume of the selected layer's 4D file
        k = sel; nv = lnv(k); newv = max(1, min(nv, lvol(k) + d));
        if nv < 2 || newv == lvol(k), return; end
        if strcmp(lfile{k}, files{1})                   % same file as the base: keep the base grid (deobliqued)
            [b, ab] = read_nifti(lfile{k}, newv, tmp); b = deoblique(b, ab); vols{k} = orient(b); nm = name_of(lfile{k});
        else
            [vols{k}, nm] = place(lfile{k}, newv);
        end
        lvol(k) = newv; lyr(k).name = lname(nm, newv, nv); refresh_all();
    end
    function select_layer(k), sel = k; refresh_all(); end
    function s = list_strings()
        s = cell(1, n);
        for k = 1:n, m = ' '; if lyr(k).vis, m = 'x'; end, s{k} = sprintf('[%s] %s  (%.2f)', m, lyr(k).name, lyr(k).alpha); end
    end
    function sl = getslice(vol, v)
        switch v
            case 1, sl = squeeze(vol(cur(1), :, :))';
            case 2, sl = squeeze(vol(:, cur(2), :))';
            otherwise, sl = squeeze(vol(:, :, cur(3)))';
        end
    end
    function refresh_all()
        xs = {(0:dim(2) - 1) * vox(2), (0:dim(1) - 1) * vox(1), (0:dim(1) - 1) * vox(1)};
        ys = {(0:dim(3) - 1) * vox(3), (0:dim(3) - 1) * vox(3), (0:dim(2) - 1) * vox(2)};
        for v = 1:3
            rgb = [];
            for k = 1:n
                if ~lyr(k).vis, continue; end
                s = getslice(vols{k}, v);
                if isempty(rgb), rgb = zeros([size(s) 3]); end
                tt = min(max((s - lyr(k).lo) / max(lyr(k).hi - lyr(k).lo, eps), 0), 1); tt(isnan(tt)) = 0;
                cm = get_cmap(lyr(k).cmap); col = reshape(cm(1 + round(tt(:) * 255), :), [size(s) 3]);
                a = lyr(k).alpha * ones(size(s)); a(isnan(s)) = 0; if lyr(k).clip, a(s <= lyr(k).lo) = 0; end
                a3 = repmat(a, [1 1 3]); rgb = rgb .* (1 - a3) + col .* a3;
            end
            if isempty(rgb), rgb = zeros(numel(ys{v}), numel(xs{v}), 3); end
            set(him(v), 'CData', rgb, 'XData', [xs{v}(1) xs{v}(end)], 'YData', [ys{v}(1) ys{v}(end)]);
            sx = {vox(2), vox(1), vox(1)}; sy = {vox(3), vox(3), vox(2)};
            set(ax(v), 'XLim', [xs{v}(1) - sx{v} / 2, xs{v}(end) + sx{v} / 2], 'YLim', [ys{v}(1) - sy{v} / 2, ys{v}(end) + sy{v} / 2]);
            px = {cur(2) - 1, cur(1) - 1, cur(1) - 1}; py = {cur(3) - 1, cur(3) - 1, cur(2) - 1};
            X = px{v} * sx{v}; Y = py{v} * sy{v};
            set(hl(1, v), 'XData', [X X], 'YData', get(ax(v), 'YLim')); set(hl(2, v), 'XData', get(ax(v), 'XLim'), 'YData', [Y Y]);
            if v == 1, set(ax(v), 'XDir', 'normal'); else, if radio, set(ax(v), 'XDir', 'reverse'); else, set(ax(v), 'XDir', 'normal'); end, end
            if v == 1, l = 'P'; r = 'A'; else, if radio, l = 'R'; r = 'L'; else, l = 'L'; r = 'R'; end, end
            if v == 3, t = 'A'; b = 'P'; else, t = 'S'; b = 'I'; end
            set(hlab(v, 1), 'String', l); set(hlab(v, 2), 'String', r); set(hlab(v, 3), 'String', t); set(hlab(v, 4), 'String', b);
        end
        set(hlist, 'String', list_strings(), 'Value', sel);
        L = lyr(sel); set(hcmap, 'Value', find(strcmp(get(hcmap, 'String'), L.cmap), 1));
        set(hlo, 'String', num2str(L.lo)); set(hhi, 'String', num2str(L.hi));
        set(halpha, 'Value', L.alpha); set(halphaval, 'String', sprintf('%.2f', L.alpha));
        set(hvis, 'Value', L.vis); set(hclip, 'Value', L.clip);
        f = P * [cur(:) - 1; 1]; fv = round(f(1:3)' + 1); w = Aras * [cur(:) - 1; 1];
        str = {sprintf('voxel (file order) [%d %d %d]   world [%.1f %.1f %.1f] mm', fv, w(1:3))};
        for k = 1:n, str{end + 1} = sprintf('%-30s %g', lyr(k).name, vols{k}(cur(1), cur(2), cur(3))); end %#ok<AGROW>
        set(hinfo, 'String', str);
    end
    function click(v)
        pt = get(ax(v), 'CurrentPoint'); x = pt(1, 1); y = pt(1, 2); act = v;
        switch v
            case 1, cur(2) = round(x / vox(2)) + 1; cur(3) = round(y / vox(3)) + 1;
            case 2, cur(1) = round(x / vox(1)) + 1; cur(3) = round(y / vox(3)) + 1;
            otherwise, cur(1) = round(x / vox(1)) + 1; cur(2) = round(y / vox(2)) + 1;
        end
        cur = max(1, min(dim, cur)); refresh_all();
    end
    function step(d)
        cur(act) = max(1, min(dim(act), cur(act) + d)); refresh_all();   % view 1 scrolls x, view 2 y, view 3 z
    end
    function key(~, e)
        switch e.Key
            case {'uparrow', 'pageup'}, step(1);
            case {'downarrow', 'pagedown'}, step(-1);
            case 'rightarrow', step(10);
            case 'leftarrow', step(-10);
            case 'g', radio = ~radio; set(hrad, 'Value', radio); refresh_all();
            case 'r', autocb();
            case 's', savepng();
            otherwise
                if strcmp(e.Character, ','), setvol(-1);
                elseif strcmp(e.Character, '.'), setvol(1);
                elseif strcmp(e.Character, '['), lyr(sel).alpha = max(0, round((lyr(sel).alpha - 0.1) * 100) / 100); refresh_all();
                elseif strcmp(e.Character, ']'), lyr(sel).alpha = min(1, round((lyr(sel).alpha + 0.1) * 100) / 100); refresh_all(); end
        end
    end
    function wheel(~, e), step(-e.VerticalScrollCount); end
    function hover(~, ~)
        fp = get(fig, 'Position'); p = get(fig, 'CurrentPoint') ./ fp(3:4);
        for v = 1:3, q = pos{v}; if p(1) >= q(1) && p(1) <= q(1) + q(3) && p(2) >= q(2) && p(2) <= q(2) + q(4), act = v; end, end
    end
    function pick(), sel = get(hlist, 'Value'); refresh_all(); end
    function addcb()
        new = nifti_pick(lastdir);
        for k = 1:numel(new), add_layer(new{k}); end
    end
    function delcb()
        if n <= 1, return; end
        vols(sel) = []; lyr(sel) = []; lfile(sel) = []; lnv(sel) = []; lvol(sel) = []; n = n - 1; sel = min(sel, n); refresh_all();
    end
    function movecb(d)
        k = sel + d; if k < 1 || k > n, return; end
        vols([sel k]) = vols([k sel]); lyr([sel k]) = lyr([k sel]); lfile([sel k]) = lfile([k sel]);
        lnv([sel k]) = lnv([k sel]); lvol([sel k]) = lvol([k sel]); sel = k; refresh_all();
    end
    function cmapcb(), items = get(hcmap, 'String'); lyr(sel).cmap = items{get(hcmap, 'Value')}; refresh_all(); end
    function rangecb()
        a = str2double(get(hlo, 'String')); b = str2double(get(hhi, 'String'));
        if ~isnan(a) && ~isnan(b) && b > a, lyr(sel).lo = a; lyr(sel).hi = b; end, refresh_all();
    end
    function alphacb(), lyr(sel).alpha = round(get(halpha, 'Value') * 100) / 100; refresh_all(); end
    function alphaeditcb()
        a = str2double(get(halphaval, 'String')); if ~isnan(a), lyr(sel).alpha = min(1, max(0, a)); end, refresh_all();
    end
    function viscb(), lyr(sel).vis = logical(get(hvis, 'Value')); refresh_all(); end
    function clipcb(), lyr(sel).clip = logical(get(hclip, 'Value')); refresh_all(); end
    function autocb(), [lyr(sel).lo, lyr(sel).hi] = auto_range(vols{sel}, sel == 1); refresh_all(); end
    function radcb(), radio = logical(get(hrad, 'Value')); refresh_all(); end
    function savepng()
        [f, d] = uiputfile('*.png', 'Save PNG', 'nifti_view.png'); if isequal(f, 0), return; end
        set(fig, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off'); print(fig, fullfile(d, f), '-dpng', '-r150');
    end
end

% ------------------------------------------------------------------------------------------------------------------
function cleanup_tmp(tmp)
if exist(tmp, 'dir'), try, rmdir(tmp, 's'); catch, end, end
end

function [lo, hi] = auto_range(v, isbase)
s = v(1:max(1, floor(numel(v) / 3e5)):end); s = s(~isnan(s) & s ~= 0);
if isempty(s), lo = 0; hi = 1; return; end
if ~isbase && min(s) >= 0 && max(s) <= 1.0001, lo = 0; hi = 1; return; end
s = sort(s); m = numel(s); lo = s(max(1, round(0.02 * m))); hi = s(max(1, round(0.998 * m)));
if hi <= lo, hi = lo + 1; end
end

function c = get_cmap(name)
t = linspace(0, 1, 256)'; z = zeros(256, 1); o = ones(256, 1);
switch name
    case 'gray', c = [t t t];
    case 'hot', c = hot(256);
    case 'jet', c = jet(256);
    case 'cool', c = cool(256);
    case 'red-yellow', c = [o t z];
    case 'blue-lightblue', c = [z t o];
    case 'green', c = [z t z];
    case 'red', c = [t z z];
    case 'blue', c = [z z t];
    otherwise, c = [t t t];
end
end

function out = resample_to(V, Av, sz, Ab)
isint = all(V(1:max(1, floor(numel(V) / 1e5)):end) == round(V(1:max(1, floor(numel(V) / 1e5)):end)));
if isint, method = 'nearest'; else, method = 'linear'; end
Fi = griddedInterpolant({single(1:size(V, 1)), single(1:size(V, 2)), single(1:size(V, 3))}, single(V), method, 'none');
out = zeros(sz, 'single'); [I, J] = ndgrid(0:sz(1) - 1, 0:sz(2) - 1); I = I(:)'; J = J(:)';
for k = 1:sz(3)
    w = Ab * [I; J; (k - 1) * ones(1, numel(I)); ones(1, numel(I))]; q = Av \ w;
    out(:, :, k) = reshape(Fi(single(q(1, :)' + 1), single(q(2, :)' + 1), single(q(3, :)' + 1)), sz(1), sz(2));
end
end

function [data, A, pixdim, nvol] = read_nifti(file, vol, tmp)
if ~exist(file, 'file'), error('File not found: %s', file); end
if endsWith(lower(file), {'.mgz', '.mgh'})
    [data, A, pixdim, nvol] = read_mgh(file, vol, tmp); return
end
if numel(file) > 3 && strcmpi(file(end - 2:end), '.gz')
    if ~exist(tmp, 'dir'), mkdir(tmp); end
    f = gunzip(file, tmp); file = f{1};
end
fid = fopen(file, 'r'); if fid < 0, error('Cannot open %s', file); end
c = onCleanup(@() fclose(fid));
fmt = 'l'; if fread(fid, 1, 'int32=>double', 0, 'l') ~= 348, fmt = 'b'; end
fseek(fid, 40, 'bof'); dim = double(fread(fid, 8, 'int16', 0, fmt));
fseek(fid, 70, 'bof'); dt = fread(fid, 1, 'int16', 0, fmt);
fseek(fid, 76, 'bof'); pixdim = double(fread(fid, 8, 'float32', 0, fmt));
fseek(fid, 108, 'bof'); voxoff = double(fread(fid, 1, 'float32', 0, fmt)); slope = double(fread(fid, 1, 'float32', 0, fmt)); inter = double(fread(fid, 1, 'float32', 0, fmt));
fseek(fid, 252, 'bof'); qcode = fread(fid, 1, 'int16', 0, fmt); scode = fread(fid, 1, 'int16', 0, fmt);
qb = double(fread(fid, 3, 'float32', 0, fmt)); qo = double(fread(fid, 3, 'float32', 0, fmt));
sx = double(fread(fid, 4, 'float32', 0, fmt)); sy = double(fread(fid, 4, 'float32', 0, fmt)); sz = double(fread(fid, 4, 'float32', 0, fmt));
types = containers.Map({2, 4, 8, 16, 64, 256, 512, 768, 1024, 1280}, {'uint8', 'int16', 'int32', 'single', 'double', 'int8', 'uint16', 'uint32', 'int64', 'uint64'});
if ~isKey(types, double(dt)), error('Unsupported NIfTI datatype %d in %s', dt, file); end
prec = types(double(dt)); nd = dim(1); d3 = dim(2:4)'; d3(d3 < 1) = 1;
bytes = struct('uint8', 1, 'int8', 1, 'int16', 2, 'uint16', 2, 'int32', 4, 'uint32', 4, 'single', 4, 'int64', 8, 'uint64', 8, 'double', 8);
nvol = max(1, dim(5) * (nd >= 4)); vol = max(1, min(vol, nvol)); fseek(fid, voxoff + (vol - 1) * prod(d3) * bytes.(prec), 'bof');
raw = fread(fid, prod(d3), [prec '=>single'], 0, fmt);
data = reshape(raw, d3); if slope ~= 0 && ~(slope == 1 && inter == 0), data = data * single(slope) + single(inter); end
if scode > 0, A = [sx'; sy'; sz'; 0 0 0 1];
elseif qcode > 0
    b = qb(1); cc = qb(2); d = qb(3); a = sqrt(max(0, 1 - b * b - cc * cc - d * d));
    R = [a * a + b * b - cc * cc - d * d, 2 * (b * cc - a * d), 2 * (b * d + a * cc); 2 * (b * cc + a * d), a * a + cc * cc - b * b - d * d, 2 * (cc * d - a * b); 2 * (b * d - a * cc), 2 * (cc * d + a * b), a * a + d * d - b * b - cc * cc];
    qf = pixdim(1); if qf == 0, qf = 1; end
    A = [R * diag([pixdim(2) pixdim(3) pixdim(4) * qf]), qo; 0 0 0 1];
else, A = diag([pixdim(2:4)', 1]);
end
pixdim = pixdim(2:4);
end

function [data, A, pixdim, nvol] = read_mgh(file, vol, tmp) % FreeSurfer .mgh / .mgz: big-endian, 284-byte header
if endsWith(lower(file), '.mgz')
    if ~exist(tmp, 'dir'), mkdir(tmp); end
    copyfile(file, fullfile(tmp, 'x.gz')); f = gunzip(fullfile(tmp, 'x.gz'), tmp); file = f{1};   % gunzip wants a .gz name
end
fid = fopen(file, 'r', 'b');
hdr = fread(fid, 7, 'int32'); dim = hdr(2:4)'; nfr = hdr(5); type = hdr(6);
fread(fid, 1, 'int16');                                      % goodRASflag
nvol = max(1, nfr); pixdim = fread(fid, 3, 'float32')'; M =reshape(fread(fid, 9, 'float32'), 3, 3); c = fread(fid, 3, 'float32');
types = {0, 'uint8'; 1, 'int32'; 3, 'single'; 4, 'int16'}; prec = types{[types{:, 1}] == type, 2};
bytes = struct('uint8', 1, 'int32', 4, 'single', 4, 'int16', 2);
fseek(fid, 284 + (min(vol, nfr) - 1) * prod(dim) * bytes.(prec), 'bof');
data = reshape(fread(fid, prod(dim), [prec '=>single']), dim);
fclose(fid);
R = M .* pixdim;                                             % direction cosines times the voxel sizes
A = [R, c - R * dim' / 2; 0 0 0 1];                          % voxel -> scanner RAS
end

function [V, A] = deoblique(V, A)                            % resample a tilted volume onto an axis-aligned grid
R = A(1:3, 1:3); s = sqrt(sum(R .^ 2));
if all(sum(abs(R ./ s) > 1e-3) == 1), return, end            % already axis-aligned (or a pure axis permutation)
sz = size(V); [i, j, k] = ndgrid([0 sz(1) - 1], [0 sz(2) - 1], [0 sz(3) - 1]);
w = A * [i(:)'; j(:)'; k(:)'; ones(1, 8)];                   % world coordinates of the 8 corners
lo = min(w(1:3, :), [], 2); n = ceil((max(w(1:3, :), [], 2) - lo)' / min(s)) + 1;
B = [min(s) * eye(3), lo; 0 0 0 1];                          % axis-aligned grid, same voxel size, covering the volume
V = resample_to(V, A, n, B); V(isnan(V)) = 0; A = B;
end
