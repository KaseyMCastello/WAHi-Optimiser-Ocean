%% iwah_geometry_drift_plots_etc.m
% Geometry-drift boxplots, optimisation diagnostics, compactness, and scaling analyses.
% X-axis uses short config IDs; a legend panel maps IDs -> names.
% Uses different colour per dimension (Δx/Δy/Δz).
%
% Requirements:
%   - RUN_DIR contains .mat files with variable H (history struct)
%   - start M: H.meta.start.M (fallback: H.iter(1).M)
%   - final M: H.meta.final.M (fallback: H.iter(end).M)
%   - CSV vs random inferred via H.meta.start.csvPath / meta.start.source
%
clear; clc;

%% ---------------- PROJECT PATHS ----------------
% Resolve repository paths from this script so execution is independent of
% the current MATLAB working directory and user account.
scriptDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(scriptDir);
runDir = fullfile(projectDir, 'data', 'runs');
figDir = fullfile(projectDir, 'fig');
resultDir = fullfile(projectDir, 'results');
addpath(fullfile(projectDir, 'src'));

if ~isfolder(figDir), mkdir(figDir); end
if ~isfolder(resultDir), mkdir(resultDir); end

saveFigures = 0;
%% ---------------- USER OPTIONS ----------------
RUN_DIR   = runDir;
FILE_GLOB = "*.mat";
AFTER_MODE = "final";  % "final" | "best"
A4_MAX_W   = 1200;

% Standard geometry map (CSV condition stem -> (configID, shortName))
% NOTE: fix spelling to match your actual stems.
STD_KEYS  = [ ...
    "4mics_Pyramid", ...
    "4mics_Tetrahedron_small", ...
    "6mics_Octahedron", ...
    "8mics_double_tetrahedorn_small", ...
    "10mics_pentagonal_antiprism", ...
    "12mics_icosahedron" ...
];
STD_NAMES = [ ...
    "Pyramid ($N_{mic}$=4)", ...
    "Tetrahedron (4)", ...
    "Octahedron (6)", ...
    "Double tetrahedron (8)", ...
    "Pentagonal antiprism (10)", ...
    "Icosahedron (12)" ...
];
STD_IDS = 1:6;

% Random config IDs: 7..11 for N=4,6,8,10,12
RAND_NS  = [4 6 8 10 12];
RAND_IDS = 7:11;
RAND_NAMES = "Random (" + string(RAND_NS) + ")";

% Visual style
STYLE = struct();
STYLE.axisFS = 10;
STYLE.titleFS = 14;

% boxplot appearance
BOX = struct();
BOX.symbol = 'o';          % outliers marker
BOX.whisker = 1.5;         % default
BOX.lineWidth = 1.2;
BOX.boxWidth  = 0.55;

% Dimension colours (different per subplot)
COL.x = [0.10 0.45 0.95];  % blue-ish
COL.y = [0.05 0.60 0.20];  % green-ish
COL.z = [0.85 0.20 0.20];  % red-ish

%% ---------------- LOAD ALL RUNS ----------------
files = dir(fullfile(RUN_DIR, FILE_GLOB));
if isempty(files)
    error('No .mat files found in "%s" matching "%s".', RUN_DIR, FILE_GLOB);
end

fprintf('Loading %d files...\n', numel(files));

R = struct([]);
for i = 1:numel(files)
    fpath = fullfile(files(i).folder, files(i).name);
    try
        S = load(fpath);
        if ~isfield(S,'H')
            warning('Skipping "%s": no H.', files(i).name);
            continue;
        end
        H = S.H;
        

        M0 = getStartM(H);
        M1_final = getFinalM(H);
        N = size(M0,1);

        % Ensure consistent N
        if size(M1_final,1) ~= N
            M1_final = M1_final(1:N,:);
        end

        [startType, condName] = inferStartConditionLocal(H);

        % AFTER mode
        switch string(AFTER_MODE)
            case "best"
                idx = pickBestIterIdx(H);
                M1 = getIterM(H, idx);
            otherwise
                M1 = M1_final;
        end

        R(end+1).file = string(files(i).name); %#ok<SAGROW>
        R(end).N = N;
        R(end).startType = string(startType);   % "csv"|"random"
        R(end).condition = string(condName);    % stem or "randomStart"
        R(end).M0 = M0;
        R(end).M1 = M1;
        R(end).H = H;   % keep full history for diagnostics
    catch ME
        warning('Failed "%s": %s', files(i).name, ME.message);
        continue;
    end
end

if isempty(R), error('No usable runs loaded.'); end
fprintf('Loaded %d runs.\n', numel(R));

%% ---------------- MAP EACH RUN -> CONFIG ID ----------------
cfgID = nan(numel(R),1);

for i = 1:numel(R)
    if R(i).startType == "csv"
        % match by condition stem
        k = find(STD_KEYS == R(i).condition, 1);
        if ~isempty(k)
            cfgID(i) = STD_IDS(k);
        else
            % unknown CSV condition -> ignore (or map to new ID if you want)
            cfgID(i) = NaN;
        end
    else
        % random: map by N to IDs 7..11
        k = find(RAND_NS == R(i).N, 1);
        if ~isempty(k)
            cfgID(i) = RAND_IDS(k);
        else
            cfgID(i) = NaN;
        end
    end
end

keep = isfinite(cfgID);
R = R(keep);
cfgID = cfgID(keep);

if isempty(R)
    error('After mapping to config IDs, no runs remained. Check STD_KEYS / RAND_NS matching.');
end

%% ---------------- BUILD DELTA VECTORS FOR BOXPLOTS ----------------
% For each run: drift per mic is (M1 - M0), then pool across mics.
dx = []; dy = []; dz = []; grp = [];

for i = 1:numel(R)
    d = R(i).M1 - R(i).M0;   % [N x 3]
    d = d(isfinite(d(:,1)) & isfinite(d(:,2)) & isfinite(d(:,3)), :);

    dx = [dx; d(:,1)]; %#ok<AGROW>
    dy = [dy; d(:,2)]; %#ok<AGROW>
    dz = [dz; d(:,3)]; %#ok<AGROW>

    grp = [grp; repmat(cfgID(i), size(d,1), 1)]; %#ok<AGROW>
end

% Categorical groups ordered 1..11
allIDs = 1:11;
g = categorical(grp, allIDs, string(allIDs), 'Ordinal', true);

%% ---------------- FIGURE: 3 BOXPLOTS + TEXTBOX LEGEND BELOW ----------------
fig = figure('Name','Geometry drift boxplots (after - before)', ...
    'Color',[1 1 1], 'Position',[80 80 980 430]);

t = tiledlayout(fig, 1, 3, 'Padding','loose', 'TileSpacing','compact');

% Reserve bottom strip for legend textbox (no extra height)
t.Units = 'normalized';
t.Position = [0.06 0.20 0.93 0.75];   % [left bottom width height] -> bottom 20% reserved



% Δx
ax1 = nexttile(t, 1);
plotBox(ax1, g, dx, COL.x, BOX, STYLE);
ylabel(ax1, '$\Delta x$ (m)', 'Interpreter','latex');
% xlabel(ax1, 'Configuration ID', 'Interpreter','latex');
title(ax1, 'X drift', 'Interpreter','latex');
ylim([-0.35 0.35])

% Δy
ax2 = nexttile(t, 2);
plotBox(ax2, g, dy, COL.y, BOX, STYLE);
ylabel(ax2, '$\Delta y$ (m)', 'Interpreter','latex');
xlabel(ax2, 'Configuration ID', 'Interpreter','latex');
title(ax2, 'Y drift', 'Interpreter','latex');
ylim([-0.35 0.35])

% Δz
ax3 = nexttile(t, 3);
plotBox(ax3, g, dz, COL.z, BOX, STYLE);
ylabel(ax3, '$\Delta z$ (m)', 'Interpreter','latex');
% xlabel(ax3, 'Configuration ID', 'Interpreter','latex');
title(ax3, 'Z drift', 'Interpreter','latex');
ylim([-0.35 0.35])

% -------- Build compact legend lines (1–2 lines) --------
stdParts  = strings(1,numel(STD_IDS));
randParts = strings(1,numel(RAND_IDS));

for k = 1:numel(STD_IDS)
    stdParts(k) = sprintf('%d=%s', STD_IDS(k), STD_NAMES(k));
end
for k = 1:numel(RAND_IDS)
    randParts(k) = sprintf('%d=%s', RAND_IDS(k), RAND_NAMES(k));
end

% Use a LaTeX-friendly separator
stdLine  = strjoin(stdParts, ';\\, ');   
randLine = strjoin(randParts, ';\\, ');

% IMPORTANT: build as char with escaped backslashes (\\textbf), not a MATLAB string literal
legendStr = sprintf('\\textbf{STD:} %s\n\\textbf{RND:} %s', stdLine, randLine);

% -------- Figure-level textbox legend --------
annotation(fig, 'textbox', [0.04 0.02 0.94 0.11], ...
    'String', legendStr, ...
    'Interpreter','latex', ...
    'FontSize', STYLE.axisFS, ...
    'EdgeColor','none', ...
    'BackgroundColor',[1 1 1], ...
    'Margin', 6, ...
    'VerticalAlignment','top', ...
    'HorizontalAlignment','center', ...
    'FitBoxToText','on', 'FontSize', 10);

% Formatting
if exist('formatLatex','file')==2
    formatLatex(ax1); formatLatex(ax2); formatLatex(ax3);
end

title(t, sprintf('Shift in mic. coordinates across optimisation runs'), ...
    'Interpreter','latex', 'FontSize', STYLE.titleFS, 'FontWeight', 'bold');

if saveFigures
    figName = 'geometry_drifts.pdf';
   exportgraphics(gcf, fullfile(figDir, figName), 'Resolution', 300) 
end

fprintf('Done.\n');

%% ===================== RUNTIME DIAGNOSTICS: WHY MAX PASS DIFFERS =====================
% Goal: compare (i) achievable max pass-rate and (ii) optimisation behaviour
% across geometries/configs: convergence speed, step sizes, and stability.

DIAG = struct();
DIAG.convFrac = 0.95;          % convergence when pass reaches 95% of its run-wise max
DIAG.useIntegerItsOnly = true; % match your trajectory logic (ignore it=0.5 shakes)
DIAG.printTable = true;

% Preallocate per-run diagnostics
nR = numel(R);
passMax    = nan(nR,1);
passFinal  = nan(nR,1);
dPass      = nan(nR,1);
nIter      = nan(nR,1);
itConv     = nan(nR,1);   % iteration index where converged (relative)
stepMean   = nan(nR,1);   % mean mic step (RMS per-iter)
stepLast   = nan(nR,1);   % last mic step
driftRMS   = nan(nR,1);   % RMS drift magnitude after-before
thr_cm     = nan(nR,1);

for i = 1:nR
    if ~isfield(R(i),'H') || ~isfield(R(i).H,'iter') || isempty(R(i).H.iter)
        continue;
    end

    H = R(i).H;

    % --- pass trajectory ---
    its  = arrayfun(@(x) x.it, H.iter(:))';
    pass = arrayfun(@(x) x.stats.passRate, H.iter(:))';  % (0..1 or 0..100?) in your batch script it was 0..1? -> check below

    % If pass is fraction, convert to percent for consistency:
    if nanmedian(pass) <= 1.5
        pass = 100 * pass;
    end

    % Optionally keep only integer iterations (ignore shake steps)
    if DIAG.useIntegerItsOnly
        isInt = abs(its - round(its)) < 1e-9;
        its  = its(isInt);
        pass = pass(isInt);
    end

    if isempty(pass)
        continue;
    end

    % --- threshold (if present) ---
    if isfield(H,'meta') && isfield(H.meta,'P') && isfield(H.meta.P,'threshold_cm')
        thr_cm(i) = H.meta.P.threshold_cm;
    end

    passMax(i)   = max(pass);
    passFinal(i) = pass(end);
    dPass(i)     = passFinal(i) - pass(1);
    nIter(i)     = numel(pass);

    % Convergence iteration: first it where pass >= convFrac * passMax
    target = DIAG.convFrac * passMax(i);
    j = find(pass >= target, 1, 'first');
    if ~isempty(j)
        itConv(i) = its(j);
    end

    % --- step-size diagnostics (how hard it had to move) ---
    % measure per-iteration mic-coordinate step RMS across mics
    % uses H.iter(k).M if available; otherwise skip.
    if isfield(H.iter(1),'M')
        Ms = cell(numel(H.iter),1);
        okM = true(numel(H.iter),1);
        for k = 1:numel(H.iter)
            if isfield(H.iter(k),'M') && ~isempty(H.iter(k).M)
                Ms{k} = H.iter(k).M;
            else
                okM(k) = false;
            end
        end

        idx = find(okM);
        if numel(idx) >= 2
            step = nan(numel(idx)-1,1);
            for kk = 2:numel(idx)
                Mprev = Ms{idx(kk-1)};
                Mcur  = Ms{idx(kk)};
                if isempty(Mprev) || isempty(Mcur), continue; end
                Nmin = min(size(Mprev,1), size(Mcur,1));
                dM = Mcur(1:Nmin,:) - Mprev(1:Nmin,:);
                step(kk-1) = sqrt(mean(dM(:).^2));  % RMS over all coords
            end
            stepMean(i) = mean(step,'omitnan');
            stepLast(i) = step(end);
        end
    end

    % --- drift magnitude (end-start) ---
    if isfield(R(i),'M0') && isfield(R(i),'M1')
        d = R(i).M1 - R(i).M0;
        driftRMS(i) = sqrt(mean(d(:).^2,'omitnan'));
    end
end

% Attach cfg IDs (already computed later, but we can compute now too)
% If cfgID exists later, you can comment this block.
cfgID_diag = nan(nR,1);
for i = 1:nR
    if R(i).startType == "csv"
        k = find(STD_KEYS == R(i).condition, 1);
        if ~isempty(k), cfgID_diag(i) = STD_IDS(k); end
    else
        k = find(RAND_NS == R(i).N, 1);
        if ~isempty(k), cfgID_diag(i) = RAND_IDS(k); end
    end
end

keepD = isfinite(cfgID_diag) & isfinite(passMax);
cfgID_diag = cfgID_diag(keepD);
passMax    = passMax(keepD);
passFinal  = passFinal(keepD);
dPass      = dPass(keepD);
nIter      = nIter(keepD);
itConv     = itConv(keepD);
stepMean   = stepMean(keepD);
stepLast   = stepLast(keepD);
driftRMS   = driftRMS(keepD);
thr_cm     = thr_cm(keepD);

% ---- Build per-config summary table ----
cfgList = unique(cfgID_diag(:))';
S = struct([]);
row = 0;

for id = cfgList
    m = (cfgID_diag == id);

    row = row + 1;
    S(row).cfgID        = id;
    S(row).nRuns        = sum(m);
    S(row).passMax_med  = median(passMax(m),'omitnan');
    S(row).passMax_iqr  = iqr(passMax(m));
    S(row).passMax_max  = max(passMax(m));
    S(row).passFinal_med= median(passFinal(m),'omitnan');
    S(row).dPass_med    = median(dPass(m),'omitnan');
    S(row).itConv_med   = median(itConv(m),'omitnan');
    S(row).stepMean_med = median(stepMean(m),'omitnan');
    S(row).stepLast_med = median(stepLast(m),'omitnan');
    S(row).driftRMS_med = median(driftRMS(m),'omitnan');

    % Name lookup for display
    if id <= 6
        S(row).name = STD_NAMES(id);
    else
        S(row).name = RAND_NAMES(id-6);
    end
end

Tdiag = struct2table(S);
Tdiag = sortrows(Tdiag, 'cfgID');

if DIAG.printTable
    disp('=== Diagnostics summary by configuration ===');
    disp(Tdiag);
end

% ---- Plot: max pass by config (box/points) ----
figD = figure('Name','Diagnostics: max/final pass and convergence','Color',[1 1 1], ...
    'Position',[120 120 950 400]);

td = tiledlayout(figD, 1, 3, 'Padding','compact', 'TileSpacing','compact');
title(td, '\textbf{Diagnostics for pass-rate discrepancies across configurations}', 'Interpreter','latex', 'FontSize', STYLE.titleFS);

gID = categorical(cfgID_diag, 1:11, string(1:11), 'Ordinal', true);

% (1) Max pass
ax = nexttile(td,1);
boxchart(ax, gID, passMax);
grid(ax,'on');
ylabel(ax,'Max pass rate (\%)','Interpreter','latex');
% xlabel(ax,'Configuration ID','Interpreter','latex');
title(ax,'Achievable max pass','Interpreter','latex');
axis(ax,'square');
box 'on';

% (2) Final pass
ax = nexttile(td,2);
boxchart(ax, gID, passFinal);
grid(ax,'on');
ylabel(ax,'Final pass rate (\%)','Interpreter','latex');
xlabel(ax,'Configuration ID','Interpreter','latex');
title(ax,'Final pass','Interpreter','latex');
axis(ax,'square');
box 'on';

% (3) Convergence iteration
ax = nexttile(td,3);
boxchart(ax, gID, itConv);
grid(ax,'on');
ylabel(ax,'Convergence iteration','Interpreter','latex');
% xlabel(ax,'Configuration ID','Interpreter','latex');
title(ax, sprintf('Convergence at %.0f\\%% of run max',100*DIAG.convFrac), 'Interpreter','latex');
axis(ax,'square');
box 'on';

if exist('formatLatex','file')==2
    formatLatex(nexttile(td,1)); formatLatex(nexttile(td,2)); formatLatex(nexttile(td,3));
end

if saveFigures
    figName = 'geometry_diagnoses.pdf';
   exportgraphics(gcf, fullfile(figDir, figName), 'Resolution', 300) 
end

% ---- Quick correlation checks (optional but useful) ----
% These help determine whether "bigger moves" or "bigger final drift" relate to higher max pass.
fprintf('\n=== Quick correlation checks (pooled runs) ===\n');
if any(isfinite(stepMean)) && any(isfinite(passMax))
    fprintf('corr(maxPass, stepMean) = %.3f\n', corr(passMax, stepMean, 'Rows','complete'));
end
if any(isfinite(driftRMS)) && any(isfinite(passMax))
    fprintf('corr(maxPass, driftRMS) = %.3f\n', corr(passMax, driftRMS, 'Rows','complete'));
end

%% ===================== GEOMETRY COMPACTNESS: BBOX SIZE BEFORE vs AFTER =====================
% Goal: quantify how "compact" the array is before vs after optimisation.
% We use axis-aligned bounding-box spans: width = x-span, height = y-span, depth = z-span.
% Summary is reported per configuration across runs.

COMPACT = struct();
COMPACT.units = "m";     % keep as metres; convert later if you prefer cm
COMPACT.printTable = true;
COMPACT.exportCSV  = true;
COMPACT.csvPath    = fullfile(resultDir, 'geometry_compactness_bbox.csv');
COMPACT.exportLatex = true;
COMPACT.latexPath   = fullfile(resultDir, 'geometry_compactness_bbox.tex');

nR = numel(R);

% Per-run spans
W0 = nan(nR,1); H0 = nan(nR,1); D0 = nan(nR,1);   % before
W1 = nan(nR,1); H1 = nan(nR,1); D1 = nan(nR,1);   % after

% Helper: bbox spans from Nx3 matrix
bboxSpans = @(M) deal( ...
    max(M(:,1)) - min(M(:,1)), ... % x-span
    max(M(:,2)) - min(M(:,2)), ... % y-span
    max(M(:,3)) - min(M(:,3))  ... % z-span
);

for i = 1:nR
    M0 = R(i).M0;
    M1 = R(i).M1;

    % Guard against NaNs or empty
    if isempty(M0) || size(M0,2)~=3 || any(~isfinite(M0(:)))
        continue;
    end
    if isempty(M1) || size(M1,2)~=3 || any(~isfinite(M1(:)))
        continue;
    end

    [W0(i), H0(i), D0(i)] = bboxSpans(M0);
    [W1(i), H1(i), D1(i)] = bboxSpans(M1);
end

% Build per-run table (useful for debugging)
TrunsCompact = table(cfgID(:), string({R.file})', [R.N]', ...
    W0, H0, D0, W1, H1, D1, ...
    'VariableNames', {'cfgID','file','Nmic','W_before','H_before','D_before','W_after','H_after','D_after'});

% Remove rows with missing values
TrunsCompact = TrunsCompact(isfinite(TrunsCompact.W_before) & isfinite(TrunsCompact.W_after), :);

% Per-config summary: median before/after
cfgList = unique(TrunsCompact.cfgID)';
S = struct([]);
row = 0;

for id = cfgList
    m = (TrunsCompact.cfgID == id);
    row = row + 1;

    % Name lookup
    if id <= 6
        cfgName = STD_NAMES(id);
    else
        cfgName = RAND_NAMES(id-6);
    end

    % Medians
    S(row).cfgID = id;
    S(row).name  = string(cfgName);
    S(row).nRuns = sum(m);

    S(row).W_before_med = median(TrunsCompact.W_before(m),'omitnan');
    S(row).H_before_med = median(TrunsCompact.H_before(m),'omitnan');
    S(row).D_before_med = median(TrunsCompact.D_before(m),'omitnan');

    S(row).W_after_med  = median(TrunsCompact.W_after(m),'omitnan');
    S(row).H_after_med  = median(TrunsCompact.H_after(m),'omitnan');
    S(row).D_after_med  = median(TrunsCompact.D_after(m),'omitnan');

    % Optional: IQR (robust spread across repeats)
    S(row).W_before_iqr = iqr(TrunsCompact.W_before(m));
    S(row).H_before_iqr = iqr(TrunsCompact.H_before(m));
    S(row).D_before_iqr = iqr(TrunsCompact.D_before(m));

    S(row).W_after_iqr  = iqr(TrunsCompact.W_after(m));
    S(row).H_after_iqr  = iqr(TrunsCompact.H_after(m));
    S(row).D_after_iqr  = iqr(TrunsCompact.D_after(m));

    % Single-number compactness: max span (largest dimension)
    S(row).maxSpan_before_med = max([S(row).W_before_med, S(row).H_before_med, S(row).D_before_med]);
    S(row).maxSpan_after_med  = max([S(row).W_after_med,  S(row).H_after_med,  S(row).D_after_med ]);

    % Ratio (<1 means became more compact, >1 means expanded)
    S(row).maxSpan_ratio = S(row).maxSpan_after_med / S(row).maxSpan_before_med;
end

Tcompact = struct2table(S);
Tcompact = sortrows(Tcompact, "cfgID");

% If you’d rather report in cm:
% cm = 100;
% varsToScale = {'W_before_med','H_before_med','D_before_med','W_after_med','H_after_med','D_after_med', ...
%                'W_before_iqr','H_before_iqr','D_before_iqr','W_after_iqr','H_after_iqr','D_after_iqr', ...
%                'maxSpan_before_med','maxSpan_after_med'};
% for v = varsToScale
%     Tcompact.(v{1}) = cm*Tcompact.(v{1});
% end
% COMPACT.units = "cm";

if COMPACT.printTable
    fprintf('\n=== Geometry compactness (axis-aligned bbox spans), units: %s ===\n', COMPACT.units);
    disp(Tcompact(:, {'cfgID','name','nRuns', ...
        'W_before_med','H_before_med','D_before_med', ...
        'W_after_med','H_after_med','D_after_med', ...
        'maxSpan_before_med','maxSpan_after_med','maxSpan_ratio'}));
end

if COMPACT.exportCSV
    writetable(Tcompact, COMPACT.csvPath);
    fprintf('Wrote: %s\n', COMPACT.csvPath);
end

if COMPACT.exportLatex
    % Simple LaTeX tabular export (no extra toolboxes)
    fid = fopen(COMPACT.latexPath,'w');
    fprintf(fid, '%% Geometry compactness table (bbox spans), units: %s\n', COMPACT.units);
    fprintf(fid, '\\begin{tabular}{r l r r r r r r r r r r}\n');
    fprintf(fid, '\\toprule\n');
    fprintf(fid, 'ID & Config & $n$ & $W_0$ & $H_0$ & $D_0$ & $W_1$ & $H_1$ & $D_1$ & $\\max_0$ & $\\max_1$ & ratio\\\\\n');
    fprintf(fid, '\\midrule\n');

    for k = 1:height(Tcompact)
        fprintf(fid, '%d & %s & %d & %.3g & %.3g & %.3g & %.3g & %.3g & %.3g & %.3g & %.3g & %.3g\\\\\n', ...
            Tcompact.cfgID(k), escapeLatexString(Tcompact.name(k)), Tcompact.nRuns(k), ...
            Tcompact.W_before_med(k), Tcompact.H_before_med(k), Tcompact.D_before_med(k), ...
            Tcompact.W_after_med(k),  Tcompact.H_after_med(k),  Tcompact.D_after_med(k), ...
            Tcompact.maxSpan_before_med(k), Tcompact.maxSpan_after_med(k), Tcompact.maxSpan_ratio(k));
    end

    fprintf(fid, '\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
    fclose(fid);
    fprintf('Wrote: %s\n', COMPACT.latexPath);
end

%% ===================== 3.5 POWER-LAW SCALING WITH MIC COUNT (CSV OUTPUT) =====================
% Outputs a CSV table (one row per metric) instead of making plots.

PLAW = struct();
PLAW.useRandomOnly = true;   % true = only random-start N sweep, false = all runs
PLAW.minN = 4;
PLAW.maxN = 12;
PLAW.bootstrap = true;
PLAW.B = 2000;               % bootstrap samples
PLAW.alpha = 0.05;           % 95% CI
PLAW.makePlots = false;      % force OFF
PLAW.exportCSV = true;
PLAW.csvPath = fullfile(resultDir, 'powerlaw_scaling_table.csv');

% Ensure output dir exists
outDir = fileparts(PLAW.csvPath);
if ~exist(outDir,'dir'); mkdir(outDir); end

% ---- Build N vector aligned to R ----
N_all = nan(numel(R),1);
for i = 1:numel(R)
    N_all(i) = R(i).N;
end

% ---- Re-derive cfg IDs for R (same mapping logic as earlier) ----
cfgID_diag_full = nan(numel(R),1);
for i = 1:numel(R)
    if R(i).startType == "csv"
        k = find(STD_KEYS == R(i).condition, 1);
        if ~isempty(k), cfgID_diag_full(i) = STD_IDS(k); end
    else
        k = find(RAND_NS == R(i).N, 1);
        if ~isempty(k), cfgID_diag_full(i) = RAND_IDS(k); end
    end
end

% ---- Align N with the already-filtered diagnostic vectors (cfgID_diag, passMax, etc.) ----
% If you kept keepD earlier, prefer using it. Otherwise we approximate with finite cfgIDs.
keepD_full = isfinite(cfgID_diag_full);
N_diag = N_all(keepD_full);

% Common keep mask across all metrics (ensures N and y are valid)
commonKeep = isfinite(N_diag) & isfinite(cfgID_diag) & isfinite(passMax);

N_use   = N_diag(commonKeep);
cfg_use = cfgID_diag(commonKeep);

% Responses
Y = struct();
Y.maxPass   = passMax(commonKeep);       % %
Y.finalPass = passFinal(commonKeep);     % %
Y.itConv    = itConv(commonKeep);        % iterations
Y.stepMean  = stepMean(commonKeep);      % m
Y.driftRMS  = driftRMS(commonKeep);      % m

% Optional random-only restriction
if PLAW.useRandomOnly
    isRandomCfg = cfg_use >= 7; % RAND IDs 7..11
else
    isRandomCfg = true(size(cfg_use));
end

% Restrict N range
isInRange = (N_use >= PLAW.minN) & (N_use <= PLAW.maxN);

maskBase = isRandomCfg & isInRange;

% Helper
fitPL = @(N,y) fitPowerLawLogLog(N,y, PLAW);

% ---- Fit all metrics and build a results table ----
keys = fieldnames(Y);

Metric   = strings(0,1);
n        = zeros(0,1);
a        = nan(0,1);
b        = nan(0,1);
a_CI_lo  = nan(0,1);
a_CI_hi  = nan(0,1);
b_CI_lo  = nan(0,1);
b_CI_hi  = nan(0,1);
R2_log   = nan(0,1);

for k = 1:numel(keys)
    nm = keys{k};
    y0 = Y.(nm);

    % Need y > 0 for log
    m = maskBase & isfinite(y0) & (y0 > 0);

    if sum(m) < 6
        warning('Power-law fit skipped for %s (too few points).', nm);
        continue;
    end

    [res, ~] = fitPL(N_use(m), y0(m));

    Metric(end+1,1) = string(nm); %#ok<SAGROW>
    n(end+1,1)      = res.n;      %#ok<SAGROW>
    a(end+1,1)      = res.a;      %#ok<SAGROW>
    b(end+1,1)      = res.b;      %#ok<SAGROW>
    R2_log(end+1,1) = res.R2;     %#ok<SAGROW>

    if isfield(res,'a_CI')
        a_CI_lo(end+1,1) = res.a_CI(1); %#ok<SAGROW>
        a_CI_hi(end+1,1) = res.a_CI(2); %#ok<SAGROW>
    else
        a_CI_lo(end+1,1) = NaN; %#ok<SAGROW>
        a_CI_hi(end+1,1) = NaN; %#ok<SAGROW>
    end

    if isfield(res,'b_CI')
        b_CI_lo(end+1,1) = res.b_CI(1); %#ok<SAGROW>
        b_CI_hi(end+1,1) = res.b_CI(2); %#ok<SAGROW>
    else
        b_CI_lo(end+1,1) = NaN; %#ok<SAGROW>
        b_CI_hi(end+1,1) = NaN; %#ok<SAGROW>
    end
end

Tpow = table(Metric, n, a, b, a_CI_lo, a_CI_hi, b_CI_lo, b_CI_hi, R2_log, ...
    'VariableNames', {'metric','n','a','b','a_CI_lo','a_CI_hi','b_CI_lo','b_CI_hi','R2_logSpace'});

disp('=== Power-law scaling summary table ===');
disp(Tpow);

if PLAW.exportCSV
    writetable(Tpow, PLAW.csvPath);
    fprintf('Wrote CSV: %s\n', PLAW.csvPath);
end

%% ---------- local helper function for power-law fit ----------
function [res, fitTbl] = fitPowerLawLogLog(N, y, PLAW)
    N = N(:); y = y(:);

    x = log(N);
    z = log(y);

    % OLS in log space: z = c0 + b*x
    X = [ones(size(x)), x];
    beta = X\z;             % [c0; b]
    c0 = beta(1);
    b  = beta(2);
    a  = exp(c0);

    % R^2 in log space
    zhat = X*beta;
    SSres = sum((z - zhat).^2);
    SStot = sum((z - mean(z)).^2);
    R2 = 1 - SSres/max(eps,SStot);

    res = struct();
    res.a = a;
    res.b = b;
    res.R2 = R2;
    res.n  = numel(N);

    fitTbl = table(N,y,x,z,zhat, 'VariableNames',{'N','y','logN','logy','logy_hat'});

    if ~PLAW.bootstrap
        return;
    end

    % Bootstrap CI for a and b
    B = PLAW.B;
    bBoot = nan(B,1);
    aBoot = nan(B,1);

    n = numel(N);
    for ii = 1:B
        idx = randi(n, n, 1);
        xb = x(idx);
        zb = z(idx);

        Xb = [ones(size(xb)), xb];
        betab = Xb\zb;
        aBoot(ii) = exp(betab(1));
        bBoot(ii) = betab(2);
    end

    res.b_CI = quantile(bBoot, [PLAW.alpha/2, 1-PLAW.alpha/2]);
    res.a_CI = quantile(aBoot, [PLAW.alpha/2, 1-PLAW.alpha/2]);
end

%% --- helper for LaTeX escaping ---
function s = escapeLatexString(str)
s = char(str);
s = strrep(s, '\', '\\textbackslash ');
s = strrep(s, '_', '\_');
s = strrep(s, '&', '\&');
s = strrep(s, '%', '\%');
s = strrep(s, '#', '\#');
s = strrep(s, '{', '\{');
s = strrep(s, '}', '\}');
s = strrep(s, '^', '\^{}');
s = strrep(s, '~', '\~{}');
end
%% ========================================================================
%% =============================== HELPERS =================================
%% ========================================================================

function plotBox(ax, g, vals, col, BOX, STYLE)
cla(ax); hold(ax,'on'); grid(ax,'on');
set(ax,'FontSize',STYLE.axisFS);

% Use boxchart for cleaner categorical control
bc = boxchart(ax, g, vals);
bc.BoxFaceColor = col;
bc.BoxFaceAlpha = 0.35;
bc.MarkerStyle  = BOX.symbol;
bc.MarkerColor  = [0 0 0];
bc.LineWidth    = BOX.lineWidth;
bc.BoxWidth     = BOX.boxWidth;

yline(ax, 0, '--k', 'Interpreter','latex');

% Keep x tick labels as IDs only
ax.XTickLabelRotation = 0;

hold(ax,'off');
axis square;
box 'on';
end

function M0 = getStartM(H)
if isfield(H,'meta') && isfield(H.meta,'start') && isfield(H.meta.start,'M')
    M0 = H.meta.start.M; return;
end
if isfield(H,'iter') && numel(H.iter)>=1 && isfield(H.iter(1),'M')
    M0 = H.iter(1).M; return;
end
error('Cannot find start mic positions (meta.start.M or iter(1).M).');
end

function M1 = getFinalM(H)
if isfield(H,'meta') && isfield(H.meta,'final') && isfield(H.meta.final,'M')
    M1 = H.meta.final.M; return;
end
if isfield(H,'iter') && numel(H.iter)>=1 && isfield(H.iter(end),'M')
    M1 = H.iter(end).M; return;
end
error('Cannot find final mic positions (meta.final.M or iter(end).M).');
end

function M = getIterM(H, idx)
if ~isfield(H,'iter') || idx<1 || idx>numel(H.iter) || ~isfield(H.iter(idx),'M')
    error('Invalid best-iter index or H.iter(idx).M missing.');
end
M = H.iter(idx).M;
end

function [startType, condName] = inferStartConditionLocal(H)
startType = "random";
condName  = "randomStart";

if isfield(H,'meta') && isfield(H.meta,'start') && isfield(H.meta.start,'csvPath')
    p = string(H.meta.start.csvPath);
    if strlength(p) > 0
        startType = "csv";
        [~,stem,~] = fileparts(p);
        condName = string(stem);
        return;
    end
end

if isfield(H,'meta') && isfield(H.meta,'start') && isfield(H.meta.start,'source')
    s = string(H.meta.start.source);
    if contains(lower(s), "csv start")
        startType = "csv";
        tok = regexp(char(s),'CSV start:\s*(.*)$','tokens','once');
        if ~isempty(tok)
            [~,stem,~] = fileparts(string(tok{1}));
            condName = string(stem);
        else
            condName = "csvStart";
        end
    end
end
end

function idx = pickBestIterIdx(H)
% Best = max passRate, tie-break by min p95Err, then max J (as in your earlier logic)
pass = arrayfun(@(x) x.stats.passRate, H.iter(:))';
p95  = arrayfun(@(x) x.stats.p95Err,    H.iter(:))';
J    = arrayfun(@(x) x.J,              H.iter(:))';

ok = isfinite(pass) & isfinite(p95) & isfinite(J);
if ~any(ok), idx = numel(H.iter); return; end

pass2 = pass; pass2(~ok) = -inf;
p952  = p95;  p952(~ok)  = inf;
J2    = J;    J2(~ok)    = -inf;

maxP = max(pass2);
cand = find(pass2==maxP);

if numel(cand)==1, idx=cand; return; end
[~,ord] = sort(p952(cand),'ascend'); cand=cand(ord);
bestP95 = p952(cand(1));
cand2 = cand(p952(cand)==bestP95);
if numel(cand2)==1, idx=cand2; return; end
[~,ord2] = sort(J2(cand2),'descend');
idx = cand2(ord2(1));
end
