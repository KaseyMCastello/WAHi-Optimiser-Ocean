%% iwah_batch_run_analyses.m
% Batch analysis + plotting for many optimiser runs saved as .mat (each containing H).
%
% Produces:
%  (1) Figure 1: Pass-% trajectories per mic-number group (4/6/8/10/12),
%      with subplots per N; TOP row = CSV standard starts, BOTTOM row = random starts.
%      Each condition gets its own colour. Median line + shaded band (IQR by default).
%
%  (2) Figure 2: “Average run” histograms of BEFORE vs AFTER error distributions per N,
%      separate CSV vs random starts.
%
%  (3) Figure 3: Accuracy vs distance-from-centre per N, separate CSV vs random.
%      Median curve + qLo–qHi band (default 10–90), BEFORE vs AFTER.
%
%  (4) Figure 4: Standard-geometry (CSV) mic scatter plots BEFORE vs AFTER
%      (overlay in one tile per geometry, dull before, bright after).
%
%  (5) Summary stats table (pooled).
%
%  (6) Delta mosaic heatmap (before→after deltas), robust colour scaling.
%
% NOTES:
%  - Each .mat contains variable H (history struct).
%  - We use integer iterations only for trajectory alignment (it = 0,1,2,...).
%    Shake records (it=0.5 etc) are ignored in trajectory plots.
%  - “Condition” for CSV starts is the CSV filename stem (without extension).
%    For random starts condition = "randomStart".
%
% REQUIREMENTS:
%  - Repository structure is preserved so the script can add src/ to the path.
%
% Ravi: edit USER OPTIONS.

clear; clc;

%% ===================== PROJECT PATHS =====================
% Resolve all repository paths from this script, not from the current
% MATLAB working directory or a user-specific absolute path.
scriptDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(scriptDir);
runDir = fullfile(projectDir, 'data', 'runs');
figDir = fullfile(projectDir, 'fig');
resultDir = fullfile(projectDir, 'results');
addpath(fullfile(projectDir, 'src'));

if ~isfolder(figDir), mkdir(figDir); end
if ~isfolder(resultDir), mkdir(resultDir); end

saveFigures = 0;
%% ===================== USER OPTIONS =====================
RUN_DIR = runDir;                      % folder containing run .mat files
FILE_GLOB = "*.mat";                   % match pattern

AFTER_MODE = "final";                  % "final" | "best" (best = max pass tie-break by p95/J)
METRIC = "ang";                        % "pos" (cm) | "ang" (deg)
IGNORE_UNUSED = true;                  % use stats.usedMask

% Thresholding for “before/after distributions”
HIST = struct();
HIST.nBins = 45;                       % histogram bin count
HIST.capMode = "4xThreshold";          % "4xThreshold" | "2xThreshold" | "none"
HIST.showCatastrophicText = true;
HIST.alphaBefore = 0.55;               % stronger pre-optimisation fill
HIST.alphaAfter  = 0.60;               % stronger post-optimisation fill
HIST.gridAlpha   = 0.10;               % keep grid subordinate to distributions
HIST.cPyrBefore = [0.35 0.35 0.35];   % grey
HIST.cPyrAfter  = [0.05 0.38 0.78];   % blue
HIST.cTetAfter  = [0.90 0.32 0.04];   % orange
HIST.cTetBefore = (1-0.42)*HIST.cTetAfter + 0.42*[1 1 1]; % light orange

% Pass-% trajectory summary band
TRAJ = struct();
TRAJ.band = "IQR";                     % "IQR" (25–75) | "10-90"
TRAJ.lineWidth = 2.0;
TRAJ.bandAlpha = 0.18;

% Accuracy-vs-distance settings
RAD = struct();
RAD.nRBins = 10;                       % radial bins for summary curves (centres)
RAD.qLo = 10;
RAD.qHi = 90;
RAD.scatterAlphaBefore = 0.25;
RAD.scatterAlphaAfter  = 0.32;
RAD.beforeLineWidth = 2.4;
RAD.afterLineWidth  = 2.4;
RAD.thresholdLineWidth = 2.0;
RAD.gridAlpha = 0.10;
RAD.capForScatter = "2xThreshold";     % "2xThreshold" | "4xThreshold" | "none"

% Figure layout
NsToShow = [4 6 8 10 12];

% ---------- Short, publication-friendly names ----------
% Map long condition keys -> short labels used EVERYWHERE (legend, row labels, tables, heatmaps)
COND_LABELS = containers.Map( ...
    [ ...
    "4mics_Pyramid", ...
    "4mics_Tetrahedron_small", ...
    "6mics_Octahedron", ...
    "8mics_double_tetrahedorn_small", ...
    "10mics_pentagonal_antiprism", ...
    "12mics_icosahedron", ...
    "randomStart" ...
    ], ...
    [ ...
    "Pyramid", ...
    "Tetrahedron", ...
    "Octahedron", ...
    "Double tetrahedron", ...
    "Pentagonal antiprism", ...
    "Icosahedron", ...
    "Random" ...
    ] ...
    );

% ---------- A4-friendly figure sizing (max width ~1200 px) ----------
FIG = struct();
FIG.maxW = 1200;
FIG.h1   = 650;    % fig1 height
FIG.h23  = 650;    % fig2/3 height
FIG.h4   = 750;    % fig4 height
FIG.hM   = 700;    % mosaic height

%% ===================== LOAD ALL RUN FILES =====================
files = dir(fullfile(RUN_DIR, FILE_GLOB));
if isempty(files)
    error('No .mat files found in RUN_DIR="%s" with pattern "%s".', RUN_DIR, FILE_GLOB);
end

Runs = struct([]);
fprintf('Loading %d files...\n', numel(files));

for i = 1:numel(files)
    fpath = fullfile(files(i).folder, files(i).name);
    try
        S = load(fpath);
        if ~isfield(S,'H')
            warning('Skipping "%s": no variable H.', files(i).name);
            continue;
        end
        H = S.H;

        if ~isfield(H,'meta') || ~isfield(H.meta,'P')
            warning('Skipping "%s": malformed H.meta.P.', files(i).name);
            continue;
        end

        P = H.meta.P;
        if ~isfield(P,'threshold_cm')
            if isfield(P,'best_cm') && isfield(P,'tol_cm')
                P.threshold_cm = P.best_cm + P.tol_cm;
            else
                error('Run "%s": P.threshold_cm missing and cannot infer from best_cm/tol_cm.', files(i).name);
            end
        end

        % Identify N
        N = size(H.meta.start.M,1);
        if isfield(H.meta,'final') && isfield(H.meta.final,'M')
            N = size(H.meta.final.M,1);
        end

        % Determine start type and geometry/condition name
        [startType, condName, csvPath] = inferStartCondition(H);

        % Choose AFTER iteration index
        iBase  = 1;
        iFinal = numel(H.iter);

        passAll = arrayfun(@(x) x.stats.passRate, H.iter(:))';
        p95All  = arrayfun(@(x) x.stats.p95Err,    H.iter(:))';
        JAll    = arrayfun(@(x) x.J,               H.iter(:))';
        iBest   = maxPassTieBreak(passAll, p95All, JAll);

        switch string(AFTER_MODE)
            case "best"
                iAfter = iBest;
            otherwise
                iAfter = iFinal;
        end

        base = H.iter(iBase);
        aft  = H.iter(iAfter);

        % Pick metric vectors
        [valsBase, labelStr] = pickMetric(base.stats, METRIC);
        [valsAft,  ~]        = pickMetric(aft.stats,  METRIC);

        % Masks
        maskBase = true(size(valsBase));
        maskAft  = true(size(valsAft));
        if IGNORE_UNUSED
            if isfield(base.stats,'usedMask'), maskBase = base.stats.usedMask(:); end
            if isfield(aft.stats,'usedMask'),  maskAft  = aft.stats.usedMask(:);  end
        end

        % Grid points (assume stable across it)
        gridPts = base.gridPts;

        % Trajectory at integer iterations only
        [itInt, passInt] = extractIntegerTrajectory(H);

        Runs(end+1).file = files(i).name; %#ok<SAGROW>
        Runs(end).path = string(fpath);
        Runs(end).N = N;
        Runs(end).startType = string(startType);   % "csv" | "random"
        Runs(end).condition = string(condName);    % long key (stem)
        Runs(end).condShort = string(condLabel(condName, COND_LABELS)); % short label for plots/tables
        Runs(end).csvPath = string(csvPath);

        Runs(end).P = P;
        Runs(end).labelStr = string(labelStr);

        Runs(end).base = base;
        Runs(end).aft  = aft;

        Runs(end).gridPts = gridPts;
        Runs(end).r = vecnorm(gridPts,2,2);

        Runs(end).valsBase = valsBase;
        Runs(end).valsAft  = valsAft;
        Runs(end).maskBase = maskBase;
        Runs(end).maskAft  = maskAft;

        Runs(end).itInt = itInt;
        Runs(end).passInt = passInt;

        % Convenience summaries
        Runs(end).passBase  = base.stats.passRate;
        Runs(end).passAfter = aft.stats.passRate;

    catch ME
        warning('Failed to load/analyse "%s": %s', files(i).name, ME.message);
        continue;
    end
end

if isempty(Runs)
    error('No usable runs were loaded.');
end

fprintf('Loaded %d runs.\n', numel(Runs));

%% ===================== GROUPING HELPERS =====================
NsAvail = unique([Runs.N]);
Ns = NsToShow(ismember(NsToShow, NsAvail));
if isempty(Ns), Ns = NsAvail; end

isCSV    = arrayfun(@(r) r.startType=="csv", Runs);
isRandom = arrayfun(@(r) r.startType=="random", Runs);

%% ===================== FIGURE 1: PASS-% TRAJECTORIES =====================
fig1 = figure('Name','Figure 1 — Pass-% trajectories by N (CSV vs Random)','Color',[1 1 1], ...
    'Position',[60 60 min(FIG.maxW,1200) FIG.h1-100]);

t1 = tiledlayout(fig1, 2, numel(Ns), 'Padding','compact', 'TileSpacing','tight');
tt = title(t1, '\textbf{Pass-rate Trajectories through Optimisation Runs}', 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter','latex');

for col = 1:numel(Ns)
    N = Ns(col);

    ax = nexttile(t1, col);
    plotPassSummary(ax, Runs, N, "csv", TRAJ); axis square;
    title(ax, sprintf('$N_{mic}=%d$', N), 'Interpreter','latex');
    cfg = configLabelForN(Runs, N, "csv");
    if strlength(cfg) > 0
        subtitle(ax, escapeLatex(cfg), 'Interpreter','latex');
    end
    ylim([10 100]); xlim([0 25]);
    ylabel('');
    xlabel(''); % no xlable for top row

    ax = nexttile(t1, numel(Ns) + col);
    plotPassSummary(ax, Runs, N, "random", TRAJ); axis square;
    ylim([0 100]); xlim([0 25]);
    title(ax, sprintf('$N_{mic}=%d$', N), 'Interpreter','latex');
    subtitle(ax, 'Random starts', 'Interpreter','latex');
    ylabel('');
    if col ~= 3
        xlabel('');
    end
end

ylabel(t1, 'Pass rate (\%)', 'Interpreter','latex');

if saveFigures && METRIC ~= "pos"
    figName = 'pass_trajectories.pdf';
   exportgraphics(gcf, fullfile(figDir, figName), 'Resolution', 300) 
   disp("pass-rate trajecotries fig saved")
end
%% ===================== FIGURE 2: “AVERAGE RUN” HISTOGRAMS BEFORE vs AFTER =====================
HIST.alphaBefore = 0.75;               % stronger pre-optimisation fill
HIST.alphaAfter  = 0.80;  

fig2 = figure('Name','Figure 2 — Average run histograms (Before vs After)','Color',[1 1 1], ...
    'Position',[80 80 min(FIG.maxW,1200) FIG.h23-100]);

t2 = tiledlayout(fig2, 2, numel(Ns), 'Padding','compact', 'TileSpacing','tight');
title(t2, sprintf('\\textbf{Average-run %s distributions: Pre vs. Post Optimisation}', ...
    escapeLatex(pickMetricName(METRIC))), 'Interpreter','latex', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter','latex');

% One xlabel only on row 2, col 3 (i.e., Ns(3) = 8 if Ns=[4 6 8 10 12])
xlCol = 3;

for col = 1:numel(Ns)
    N = Ns(col);

    % -------- TOP ROW: CSV starts --------
    ax = nexttile(t2, col);
    plotAverageHistogram(ax, Runs, N, "csv", METRIC, HIST);

    title(ax, sprintf('$N_{mic}=%d$', N), 'Interpreter','latex');

    cfg = configLabelForN(Runs, N, "csv");
    if strlength(cfg) > 0
        subtitle(ax, escapeLatex(cfg), 'Interpreter','latex');
    end

    if METRIC == "pos"
        ylim([0 0.25]);
    else
        % ylim([0 0.5]);
    end
    if METRIC == "ang"
        xlim([0 15]);
        % ylim([0 0.8])
    end

    ylabel(ax, '');
    xlabel(ax,'');
    axis(ax,'square');

    % -------- BOTTOM ROW: Random starts --------
    ax = nexttile(t2, numel(Ns)+col);
    plotAverageHistogram(ax, Runs, N, "random", METRIC, HIST);
    if METRIC == "pos"
        ylim([0 0.25]);
    else
        % ylim([0 0.5]);
    end
    title(ax, sprintf('$N_{mic}=%d$', N), 'Interpreter','latex');
    subtitle(ax, 'Random starts', 'Interpreter','latex');

    ylabel(ax, '');

    if col == xlCol
        xlabel(ax, pickMetricName(METRIC), 'Interpreter','latex');
    else
        xlabel(ax,'');
    end

    if METRIC == "ang"
        xlim([0 15]);
        % ylim([0 0.8])
    end
    axis(ax,'square');
end

ylabel(t2, 'Mean probability (across runs)', 'Interpreter','latex');

if saveFigures
    figName = strcat(['histograms_' char(METRIC) '_error.pdf']);
   exportgraphics(gcf, fullfile(figDir, figName), 'Resolution', 300) 
end

%% ===================== FIGURE 3: ACCURACY vs DISTANCE FROM CENTRE =====================
HIST.alphaBefore = 0.55;               % stronger pre-optimisation fill
HIST.alphaAfter  = 0.60;  

fig3 = figure('Name','Figure 3 — Accuracy vs distance from centre', ...
    'Color',[1 1 1], 'Position',[90 90 min(FIG.maxW,1200) FIG.h23-100]);

t3 = tiledlayout(fig3, 2, numel(Ns), ...
    'Padding','compact', 'TileSpacing','compact');

title(t3, sprintf('\\textbf{%s vs distance from array centre (median + %d--%d\\%% band)}', ...
    escapeLatex(pickMetricName(METRIC)), RAD.qLo, RAD.qHi), ...
    'Interpreter','latex', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter','latex');

% Column where we place the single x-label (middle column)
xlCol = ceil(numel(Ns)/2);

for col = 1:numel(Ns)
    N = Ns(col);

    % -------- TOP ROW: CSV --------
    ax = nexttile(t3, col);
    plotRadialAccuracy(ax, Runs, N, "csv", METRIC, RAD, HIST);
    title(ax, sprintf('$N_{mic}=%d$', N), 'Interpreter','latex');
    cfg = configLabelForN(Runs, N, "csv");
    if strlength(cfg) > 0
        subtitle(ax, escapeLatex(cfg), 'Interpreter','latex');
    end

    axis(ax,'square');
    if METRIC == "pos"
        % ylim([0 35]);
    elseif METRIC == "ang"
        % ylim([0 10]);
    end

    ylabel(ax,'');

    xlabel(ax,'');  % no x-label in top row

    % -------- BOTTOM ROW: RANDOM --------
    ax = nexttile(t3, numel(Ns) + col);
    plotRadialAccuracy(ax, Runs, N, "random", METRIC, RAD, HIST);
    title(ax, sprintf('$N_{mic}=%d$', N), 'Interpreter','latex');
    subtitle(ax, 'Random starts', 'Interpreter','latex');
    axis(ax,'square');
    if METRIC == "pos"
        % ylim([0 35]);
    elseif METRIC == "ang"
        % ylim([0 10]);
    end

    ylabel(ax,'');

    % Single shared x-label
    if col == xlCol
        xlabel(ax, 'Distance from array centre (m)', 'Interpreter','latex');
    else
        xlabel(ax,'');
    end
end

ylabel(t3, pickMetricName(METRIC), 'Interpreter','latex');

if saveFigures
    figName = strcat(['accuracy_' char(METRIC) '_distance.pdf']);
   exportgraphics(gcf, fullfile(figDir, figName), 'Resolution', 300) 
end

%% ========================================================================
%% Single-geometry BEFORE vs AFTER (Figure-4 style, but only one geometry)
%% ========================================================================

% -------------------- USER SELECT --------------------
% Use either a SHORT label (from COND_LABELS map), e.g.:
%   "Icosahedron", "Octahedron", "Pentagonal antiprism", "Pyramid", ...
% OR use the RAW condition key stem (from CSV filename), e.g.:
%   "12mics_icosahedron", "6mics_Octahedron", ...

SELECT_GEOM = "Pyramid";


% -------------------- FILTER CSV RUNS --------------------
isCSV = arrayfun(@(r) r.startType=="csv", Runs);
csvRuns = Runs(isCSV);

if isempty(csvRuns)
    error('No CSV runs found in Runs.');
end

% Build available lists (both raw and short)
rawConds   = unique([csvRuns.condition]);
shortConds = unique([csvRuns.condShort]);

% Resolve SELECT_GEOM to a raw condition key (csvRuns.condition)
select = string(SELECT_GEOM);

% Try: exact match to raw condition
hitRaw = rawConds(rawConds == select);

% Try: exact match to short label
hitShort = rawConds(arrayfun(@(c) any([csvRuns.condition]==c & [csvRuns.condShort]==select), rawConds));

% Robust alternative: search within csvRuns directly
if isempty(hitRaw) && isempty(hitShort)
    % look for rows where either condition OR condShort matches
    idx = find([csvRuns.condition]==select | [csvRuns.condShort]==select, 1);
    if isempty(idx)
        fprintf('Available SHORT geometries:\n');
        disp(shortConds(:));
        fprintf('Available RAW condition keys:\n');
        disp(rawConds(:));
        error('Could not find geometry "%s" in Runs.', select);
    end
    condKey = csvRuns(idx).condition;
elseif ~isempty(hitRaw)
    condKey = hitRaw(1);
else
    condKey = hitShort(1);
end

% Pick representative run for that geometry
run = pickRepresentativeRun(csvRuns, condKey);

% Convenient labels
condShort = condLabel(condKey, COND_LABELS);
N = run.N;

% -------------------- STYLE / THRESHOLDS --------------------
P        = run.P;
thr      = P.threshold_cm;
cmax     = 2*thr;
catThr   = 4*thr;

gridPts  = run.gridPts;
labelStr = run.labelStr;

% Before/After values + masks
valsBase = run.valsBase;  maskBase = run.maskBase;
valsAft  = run.valsAft;   maskAft  = run.maskAft;

% -------------------- FIGURE 4 --------------------
fig = figure('Name', sprintf('Single geometry: %s (N=%d)', condShort, N), ...
    'Color',[1 1 1], 'Position',[200 200 900 400]);

t = tiledlayout(fig, 1, 2, 'Padding','compact', 'TileSpacing','compact');

% title(t, sprintf('\\textbf{%s} --- $N_{mic}=%d$', ...
%     escapeLatex(condShort), N), 'Interpreter','latex');

% -------- BEFORE --------
ax1 = nexttile(t, 1);
plot3D(ax1, P, run.base.M, gridPts, valsBase, maskBase, labelStr, thr, cmax, catThr);
title(ax1, sprintf('Start Config: pass=%.1f\\%%', 100*run.passBase), 'Interpreter','latex');
xlim(ax1, [-0.5 2]); ylim(ax1, [-2 2]); zlim(ax1, [-0.5 2]);
colorbar off;
formatLatex(gca);
grid minor
% -------- AFTER --------
ax2 = nexttile(t, 2);
plot3D(ax2, P, run.aft.M, gridPts, valsAft, maskAft, labelStr, thr, cmax, catThr);
title(ax2, sprintf('Post-optimisation: pass=%.1f\\%%', 100*run.passAfter), 'Interpreter','latex');
xlim(ax2, [-0.5 2]); ylim(ax2, [-2 2]); zlim(ax2, [-0.5 2]);
formatLatex(gca);
grid minor
% Optional: consistent viewpoint
view(ax1, 3);
view(ax2, 3);

% Optional: add a small subtitle with thresholds (safe interpreter)
sgtitle('\textbf{Optimisation Run - Demonstration}', 'FontSize', 14, ...
    'FontWeight','bold', 'Interpreter', 'latex');

if saveFigures && METRIC == "pos"
    figName = 'optimisation_demo.pdf';
   exportgraphics(gcf, fullfile(figDir, figName), 'Resolution', 300) 
end
%% ===================== SUMMARY STATS TABLE (POOLED) =====================
thr0 = median(arrayfun(@(r) r.P.threshold_cm, Runs), 'omitnan'); %#ok<NASGU>
poolVectors = @(R) pooledBeforeAfterVectors(R, METRIC);

rows = struct([]);
row = 0;

csvRuns = Runs(arrayfun(@(r) r.startType=="csv", Runs));
if ~isempty(csvRuns)
    Ns_csv = unique([csvRuns.N]);
    for iN = 1:numel(Ns_csv)
        N = Ns_csv(iN);
        Rn = csvRuns([csvRuns.N] == N);
        conds = unique([Rn.condition]);

        for ic = 1:numel(conds)
            cond = conds(ic);
            Rc = Rn([Rn.condition] == cond);
            if isempty(Rc), continue; end

            [vB, vA, thr] = poolVectors(Rc);

            row = row + 1;
            rows(row).startType = "csv";
            rows(row).condition = condLabel(cond, COND_LABELS);  % <- SHORT
            rows(row).N = N;
            rows(row).thr = thr;

            rows(row).pass_before   = 100 * mean(vB <= thr, 'omitnan');
            rows(row).mean_before   = mean(vB, 'omitnan');
            rows(row).median_before = median(vB, 'omitnan');
            qB = prctile(vB, [25 75]); rows(row).iqr_before = qB(2) - qB(1);
            rows(row).p95_before    = prctile(vB, 95);

            rows(row).pass_after    = 100 * mean(vA <= thr, 'omitnan');
            rows(row).mean_after    = mean(vA, 'omitnan');
            rows(row).median_after  = median(vA, 'omitnan');
            qA = prctile(vA, [25 75]); rows(row).iqr_after = qA(2) - qA(1);
            rows(row).p95_after     = prctile(vA, 95);
        end
    end
end

rndRuns = Runs(arrayfun(@(r) r.startType=="random", Runs));
if ~isempty(rndRuns)
    Ns_rnd = unique([rndRuns.N]);
    for iN = 1:numel(Ns_rnd)
        N = Ns_rnd(iN);
        Rn = rndRuns([rndRuns.N] == N);
        if isempty(Rn), continue; end

        [vB, vA, thr] = poolVectors(Rn);

        row = row + 1;
        rows(row).startType = "random";
        rows(row).condition = condLabel("randomStart", COND_LABELS); % <- SHORT ("Random")
        rows(row).N = N;
        rows(row).thr = thr;

        rows(row).pass_before   = 100 * mean(vB <= thr, 'omitnan');
        rows(row).mean_before   = mean(vB, 'omitnan');
        rows(row).median_before = median(vB, 'omitnan');
        qB = prctile(vB, [25 75]); rows(row).iqr_before = qB(2) - qB(1);
        rows(row).p95_before    = prctile(vB, 95);

        rows(row).pass_after    = 100 * mean(vA <= thr, 'omitnan');
        rows(row).mean_after    = mean(vA, 'omitnan');
        rows(row).median_after  = median(vA, 'omitnan');
        qA = prctile(vA, [25 75]); rows(row).iqr_after = qA(2) - qA(1);
        rows(row).p95_after     = prctile(vA, 95);
    end
end

T = struct2table(rows);
T.startType = string(T.startType);
T.condition = string(T.condition);
T = sortrows(T, {'startType','N','condition'});
disp(T);
saveCSV = 1;
if saveCSV
    writetable(T, fullfile(resultDir, ['summary_stats_' char(METRIC) '.csv']))
end
%% ===================== DELTA MOSAIC (BEFORE→AFTER) =====================
% Robust colour scaling (95th percentile of |Δ|), symmetric about 0.
% Uses SHORT condition names everywhere.
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
rwo_ids = [STD_IDS RAND_IDS];


MOSAIC = struct();
MOSAIC.capMult    = 4;
MOSAIC.showTail   = true;
MOSAIC.fontSize   = 12;
MOSAIC.robustPct  = 95;
MOSAIC.pdfSizeIn  = [6 4.5]; % fixed page size keeps POS and ANG PDFs identical

useCap = (string(METRIC) == "pos");

D = makeDeltaTableFromRuns(Runs, METRIC, useCap, MOSAIC.capMult, COND_LABELS);

colNames = ["dPass_pp","dMean","dMedian","dIQR","dP95"];
niceCol  = ["$\Delta$ pass","$\Delta$ mean","$\Delta$ median","$\Delta$ IQR","$\Delta$ p95"];
if MOSAIC.showTail
    colNames = [colNames, "dTail_pp"];
    niceCol  = [niceCol, "$\Delta$ tail"];
end

Z = zeros(height(D), numel(colNames));
for j = 1:numel(colNames)
    v = D.(colNames(j));
    if colNames(j) == "dPass_pp"
        Z(:,j) = v;
    else
        Z(:,j) = -v; % improvement => positive
    end
end

rowLbl = D.rowLabel;

absZ = abs(Z(:));
absZ = absZ(isfinite(absZ) & absZ > 0);
if isempty(absZ)
    clim = 1;
else
    clim = prctile(absZ, MOSAIC.robustPct);
    if ~isfinite(clim) || clim <= 0, clim = max(absZ); end
end

figM = figure('Name','Delta mosaic (improvement map)','Color',[1 1 1], ...
    'Position',[150 150 600 450]);

imagesc(Z);
axis tight; axis ij;

set(gca,'XTick',1:numel(niceCol), ...
    'XTickLabel',niceCol, ...
    'YTick',1:numel(rwo_ids), ...
    'YTickLabel',rwo_ids, ...
    'TickLabelInterpreter','latex', ...
    'FontSize',MOSAIC.fontSize);

ylabel('Configuration ID');
colormap(divergingBlueRedGreen(256));
caxis([-clim, +clim]);
c = colorbar;
c.Label.String = sprintf('improvement');
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
% tt = title(sprintf('\Delta (after vs before): green=improvement, red=worsening, blue=neutral \;|\; METRIC=%s', ...
%     escapeLatex(string(METRIC))), 'Interpreter','latex');

if useCap
    st = subtitle(sprintf('POS deltas', MOSAIC.capMult), ...
        'Interpreter','latex');
else
    st = subtitle('ANG deltas', 'Interpreter','latex');
end

for i = 1:size(Z,1)
    for j = 1:size(Z,2)
        if colNames(j) == "dPass_pp"
            txt = sprintf('%+.1f', D.dPass_pp(i));
        elseif colNames(j) == "dTail_pp"
            txt = sprintf('%+.1f', D.dTail_pp(i));
        else
            txt = sprintf('%+.2g', D.(colNames(j))(i));
        end
        text(j, i, txt, 'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'Color','w','FontSize',MOSAIC.fontSize-1,'FontWeight','bold', 'Interpreter','none');
    end
end
formatLatex(gca);
grid off
saveFigures = 1;
if saveFigures
    figName = strcat([char(METRIC) '_improvements.pdf']);
    set(figM, ...
        'PaperUnits','inches', ...
        'PaperSize',MOSAIC.pdfSizeIn, ...
        'PaperPosition',[0 0 MOSAIC.pdfSizeIn], ...
        'PaperPositionMode','manual');
    print(figM, fullfile(figDir, figName), '-dpdf', '-painters');
end
saveFigures = 0;
%% ========================================================================
%% ============================ LOCAL FUNCTIONS ===========================
%% ========================================================================
function g = configLabelForN(Runs, N, startType)
% Returns a compact config label for the given N and startType.
sel = arrayfun(@(r) r.N==N && r.startType==startType, Runs);
R = Runs(sel);

if isempty(R)
    g = "";
    return;
end

if string(startType) == "csv"
    geoms = unique(string({R.condShort}), 'stable');
    geoms = geoms(geoms ~= "");
    if isempty(geoms)
        g = "";
    else
        g = strjoin(cellstr(geoms), " / ");
    end
else
    g = "Random";
end
end

function s = condLabel(condName, COND_LABELS)
% Map long condition key -> short label (fallback to original)
condName = string(condName);
if isKey(COND_LABELS, condName)
    s = string(COND_LABELS(condName));
else
    s = condName;
end
end

function out = escapeLatex(s)
% Minimal LaTeX escaping for titles (not for tick labels you want verbatim)
s = string(s);
out = replace(s, "\", "\\");
out = replace(out, "_", "\_");
out = replace(out, "%", "\%");
out = replace(out, "&", "\&");
out = replace(out, "#", "\#");
out = replace(out, "{", "\{");
out = replace(out, "}", "\}");
out = replace(out, "^", "\^{}");
out = replace(out, "~", "\~{}");
out = char(out);
end

function [startType, condName, csvPath] = inferStartCondition(H)
csvPath = "";
condName = "";
startType = "random";

if isfield(H,'meta') && isfield(H.meta,'start') && isfield(H.meta.start,'csvPath')
    p = string(H.meta.start.csvPath);
    if strlength(p) > 0
        startType = "csv";
        csvPath = p;
        [~,stem,~] = fileparts(p);
        condName = stem;
        return;
    end
end

if isfield(H,'meta') && isfield(H.meta,'start') && isfield(H.meta.start,'source')
    s = string(H.meta.start.source);
    if contains(lower(s),"csv start")
        startType = "csv";
        tok = regexp(char(s),'CSV start:\s*(.*)$','tokens','once');
        if ~isempty(tok)
            csvPath = string(tok{1});
            [~,stem,~] = fileparts(csvPath);
            condName = string(stem);
        else
            condName = "csvStart";
        end
        return;
    end
end

startType = "random";
condName = "randomStart";
end

function [itInt, passInt] = extractIntegerTrajectory(H)
itAll = arrayfun(@(x) x.it, H.iter(:))';
passAll = arrayfun(@(x) x.stats.passRate, H.iter(:))';

isInt = abs(itAll - round(itAll)) < 1e-9;
itInt = itAll(isInt);
passInt = passAll(isInt);

[itInt, idx] = sort(itInt);
passInt = passInt(idx);

if isempty(itInt) || itInt(1) ~= 0
    itInt = [0; itInt(:)];
    passInt = [passAll(1); passInt(:)];
end
end

function plotPassSummary(ax, Runs, N, startType, TRAJ)
grid(ax,'on'); hold(ax,'on');

sel = arrayfun(@(r) r.N==N && r.startType==startType, Runs);
R = Runs(sel);

if isempty(R)
    text(ax,0.5,0.5,'No runs','Units','normalized','HorizontalAlignment','center');
    axis(ax,'off');
    return;
end

% Use SHORT labels for grouping. Keep the two four-microphone canonical
% geometries in a fixed order so their colours and annotations cannot
% change with file-loading order.
conds = unique([R.condShort], 'stable');
isFourMicCanonical = N == 4 && string(startType) == "csv";
if isFourMicCanonical
    preferredOrder = ["Pyramid", "Tetrahedron"];
    preferredOrder = preferredOrder(ismember(preferredOrder, conds));
    conds = [preferredOrder, conds(~ismember(conds, preferredOrder))];
end
cols  = lines(max(1,numel(conds)));

maxIt = max(arrayfun(@(r) max(r.itInt), R));
x = (0:maxIt)';

switch string(TRAJ.band)
    case "10-90"
        qlo = 10; qhi = 90;
    otherwise
        qlo = 25; qhi = 75;
end

hLine = gobjects(numel(conds),1);  % store ONLY line handles for legend

for c = 1:numel(conds)
    condS = conds(c);
    Rc = R([R.condShort]==condS);
    if isempty(Rc), continue; end

    Y = nan(numel(x), numel(Rc));
    for k = 1:numel(Rc)
        y = alignTrajectoryToX(Rc(k).itInt, 100*Rc(k).passInt, x);
        Y(:,k) = y;
    end

    med = median(Y,2,'omitnan');
    lo  = prctile(Y,qlo,2);
    hi  = prctile(Y,qhi,2);

    % Band (exclude from legend)
    hBand = fill(ax, [x; flipud(x)], [lo; flipud(hi)], cols(c,:), ...
        'FaceAlpha', TRAJ.bandAlpha, 'EdgeColor','none');
    hBand.HandleVisibility = 'off';  % <-- key: legend ignores this
    

    % Line (this is what legend should label)
    hLine(c) = plot(ax, x, med, 'LineWidth', TRAJ.lineWidth, 'Color', cols(c,:));

end

% Explicit in-panel colour key requested during peer review. Text is placed
% in normalised axes coordinates so it remains stable if limits change.
if isFourMicCanonical
    pyramidIdx = find(conds == "Pyramid", 1);
    tetrahedronIdx = find(conds == "Tetrahedron", 1);

    if ~isempty(pyramidIdx)
        text(ax, 0.57, 0.18, '-- Pyramid', ...
            'Units', 'normalized', ...
            'Color', cols(pyramidIdx, :), ...
            'Interpreter', 'none', ...
            'FontWeight', 'bold', ...
            'FontSize', 9, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle');
    end
    if ~isempty(tetrahedronIdx)
        text(ax, 0.57, 0.08, '-- Tetrahedron', ...
            'Units', 'normalized', ...
            'Color', cols(tetrahedronIdx, :), ...
            'Interpreter', 'none', ...
            'FontWeight', 'bold', ...
            'FontSize', 9, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle');
    end
end

xlabel(ax,'Iteration', 'Interpreter','latex');
ylabel(ax,'Pass rate (\%)','Interpreter','latex');
ylim(ax,[0 100]);

% Legend: only lines + matching labels
ok = isgraphics(hLine);
% legend(ax, hLine(ok), cellstr(conds(ok)), 'Interpreter','latex', 'Location','south');

formatLatex(ax);
hold(ax,'off');
end

function yAligned = alignTrajectoryToX(it, y, x)
yAligned = nan(size(x));
it = it(:); y = y(:);

for i = 1:numel(it)
    idx = find(x==it(i), 1);
    if ~isempty(idx)
        yAligned(idx) = y(i);
    end
end

last = nan;
for i = 1:numel(yAligned)
    if isfinite(yAligned(i))
        last = yAligned(i);
    else
        if isfinite(last)
            yAligned(i) = last;
        end
    end
end
end

function plotAverageHistogram(ax, Runs, N, startType, METRIC, HIST)
% Filled "average-run" histograms (mean probability across runs).
% Special case: N=4 & CSV -> Pyramid (grey/blue) + Tetrahedron (light-orange/orange).

HIST = ensureHistogramStyle(HIST);

grid(ax,'on'); hold(ax,'on');
ax.GridAlpha = HIST.gridAlpha;
ax.MinorGridAlpha = 0.5 * HIST.gridAlpha;

sel = arrayfun(@(r) r.N==N && r.startType==startType, Runs);
R = Runs(sel);

if isempty(R)
    text(ax,0.5,0.5,'No runs','Units','normalized','HorizontalAlignment','center');
    axis(ax,'off'); return;
end

thr = median(arrayfun(@(r) r.P.threshold_cm, R), 'omitnan');

% Cap logic (NO subtitle printed)
switch string(HIST.capMode)
    case "2xThreshold", cap = 2*thr;
    case "4xThreshold", cap = 4*thr;
    otherwise,          cap = inf;
end

% ----- Build common edges (shared across everything drawn in this panel) -----
allVals = [];
for k = 1:numel(R)
    [vB, vA] = getRunMetricVectors(R(k), METRIC);
    allVals = [allVals; vB(:); vA(:)]; %#ok<AGROW>
end
allVals = allVals(isfinite(allVals));
if isfinite(cap), allVals = allVals(allVals<=cap); end

if isempty(allVals)
    text(ax,0.5,0.5,'No finite values','Units','normalized','HorizontalAlignment','center');
    axis(ax,'off'); return;
end

vMax  = max(allVals);
edges = linspace(0, vMax, HIST.nBins+1);
centres = (edges(1:end-1)+edges(2:end))/2;

% ----- Colours -----
cPyr_before = HIST.cPyrBefore;
cPyr_after  = HIST.cPyrAfter;
cTet_before = HIST.cTetBefore;
cTet_after  = HIST.cTetAfter;

alphaBefore = HIST.alphaBefore;
alphaAfter  = HIST.alphaAfter;

% ----- Draw -----
if string(startType)=="csv" && N==4
    % Split N=4 CSV into Pyramid vs Tetrahedron
    shorts = string({R.condShort});
    conds  = string({R.condition});

    isPyr = (shorts=="Pyramid")     | contains(lower(conds), "pyramid");
    isTet = (shorts=="Tetrahedron") | contains(lower(conds), "tetra");

    Rp = R(isPyr);
    Rt = R(isTet);

    % Fall back safely if split fails
    if isempty(Rp) || isempty(Rt)
        [pB, pA] = pooledHistProb(R, METRIC, edges, cap);
        drawFilled(ax, centres, pB, cPyr_before, alphaBefore);
        drawFilled(ax, centres, pA, cPyr_after,  alphaAfter);
    else
        [pB_p, pA_p] = pooledHistProb(Rp, METRIC, edges, cap);
        [pB_t, pA_t] = pooledHistProb(Rt, METRIC, edges, cap);

        % Order: before then after, so after sits on top
        drawFilled(ax, centres, pB_p, cPyr_before, alphaBefore);
        drawFilled(ax, centres, pA_p, cPyr_after,  alphaAfter);

        drawFilled(ax, centres, pB_t, cTet_before, alphaBefore);
        drawFilled(ax, centres, pA_t, cTet_after,  alphaAfter);
    end
else
    % Default: pooled before/after (grey/blue)
    [pB, pA] = pooledHistProb(R, METRIC, edges, cap);
    drawFilled(ax, centres, pB, cPyr_before, alphaBefore);
    drawFilled(ax, centres, pA, cPyr_after,  alphaAfter);
end

% Threshold line
xline(ax, thr, '--k', 'thr');

% Labels (you will blank these at the layout level as needed)
xlabel(ax, pickMetricName(METRIC), 'Interpreter','latex');
ylabel(ax, 'Mean probability (across runs)', 'Interpreter','latex');

% Catastrophic textbox (optional)
if isfield(HIST,'showCatastrophicText') && HIST.showCatastrophicText && isfinite(cap)
    nCatB = 0; nTotB = 0;
    nCatA = 0; nTotA = 0;
    for k = 1:numel(R)
        [vB, vA] = getRunMetricVectors(R(k), METRIC);
        nTotB = nTotB + numel(vB);
        nTotA = nTotA + numel(vA);
        nCatB = nCatB + sum(vB>cap);
        nCatA = nCatA + sum(vA>cap);
    end

    fracB = nCatB / max(1,nTotB);
    fracA = nCatA / max(1,nTotA);

    % txt = sprintf(['Catastrophic ($>$ cap)\n' ...
    %                'before: %.1f\\%% (%d/%d)\n' ...
    %                'after:  %.1f\\%% (%d/%d)'], ...
    %                100*fracB, nCatB, nTotB, ...
    %                100*fracA, nCatA, nTotA);
    %
    % text(ax, 0.98, 0.98, txt, ...
    %     'Units','normalized', ...
    %     'HorizontalAlignment','right', ...
    %     'VerticalAlignment','top', ...
    %     'FontSize',9, ...
    %     'Interpreter','latex', ...
    %     'BackgroundColor',[1 1 1 0.85], ...
    %     'EdgeColor',[0.75 0.75 0.75]);
end

formatLatex(ax);
hold(ax,'off');
end

% ---- helpers (keep as local functions in the same file) ----
function drawFilled(ax, centres, p, col, alpha)
p = p(:)'; centres = centres(:)';
fill(ax, [centres, fliplr(centres)], [p, zeros(size(p))], col, ...
    'FaceAlpha', alpha, 'EdgeColor','none');
end

function [pB, pA] = pooledHistProb(R, METRIC, edges, cap)
HB = nan(numel(edges)-1, numel(R));
HA = nan(numel(edges)-1, numel(R));
for k = 1:numel(R)
    [vB, vA] = getRunMetricVectors(R(k), METRIC);
    if isfinite(cap)
        vB = vB(vB<=cap);
        vA = vA(vA<=cap);
    end
    HB(:,k) = histcounts(vB, edges, 'Normalization','probability')';
    HA(:,k) = histcounts(vA, edges, 'Normalization','probability')';
end
pB = mean(HB,2,'omitnan');
pA = mean(HA,2,'omitnan');
end


function [vB, vA] = getRunMetricVectors(run, METRIC)
[vB0, ~] = pickMetric(run.base.stats, METRIC);
[vA0, ~] = pickMetric(run.aft.stats,  METRIC);

mB = run.maskBase(:);
mA = run.maskAft(:);

vB = vB0(mB);
vA = vA0(mA);

vB = vB(isfinite(vB));
vA = vA(isfinite(vA));
end

function idx = maxPassTieBreak(pass, p95, J)
pass = pass(:);
p95  = p95(:);
J    = J(:);

ok = isfinite(pass) & isfinite(p95) & isfinite(J);
if ~any(ok)
    idx = numel(pass);
    return;
end

pass2 = pass; pass2(~ok) = -inf;
p952  = p95;  p952(~ok)  = inf;
J2    = J;    J2(~ok)    = -inf;

maxP = max(pass2);
cand = find(pass2 == maxP);

if numel(cand) == 1
    idx = cand;
    return;
end

[~, ord] = sort(p952(cand), 'ascend');
cand = cand(ord);
bestP95 = p952(cand(1));
cand2 = cand(p952(cand) == bestP95);

if numel(cand2) == 1
    idx = cand2;
    return;
end

[~, ord2] = sort(J2(cand2), 'descend');
idx = cand2(ord2(1));
end

function [vals, labelStr] = pickMetric(stats, METRIC)
m = string(METRIC);

if m == "ang"
    if ~isfield(stats,'ang_deg')
        error('pickMetric:MissingField','stats.ang_deg not found.');
    end
    vals = stats.ang_deg(:);
    labelStr = "Angular error (deg)";
else
    if ~isfield(stats,'err_cm')
        error('pickMetric:MissingField','stats.err_cm not found.');
    end
    vals = stats.err_cm(:);
    labelStr = "Position error (cm)";
end
end

function nameStr = pickMetricName(METRIC)
m = string(METRIC);
if m == "ang"
    nameStr = "Angular error (deg)";
else
    nameStr = "Position error (cm)";
end
end

function plotRadialAccuracy(ax, Runs, N, startType, METRIC, RAD, HIST)
HIST = ensureHistogramStyle(HIST);

grid(ax,'on'); hold(ax,'on');
ax.GridAlpha = RAD.gridAlpha;
ax.MinorGridAlpha = 0.5 * RAD.gridAlpha;

sel = arrayfun(@(r) r.N==N && r.startType==startType, Runs);
R = Runs(sel);

if isempty(R)
    text(ax,0.5,0.5,'No runs','Units','normalized','HorizontalAlignment','center');
    axis(ax,'off');
    return;
end

thr = median(arrayfun(@(r) r.P.threshold_cm, R), 'omitnan');

cap = inf;
switch string(RAD.capForScatter)
    case "2xThreshold", cap = 2*thr;
    case "4xThreshold", cap = 4*thr;
    otherwise, cap = inf;
end

if string(startType) == "csv"
    conds = unique([R.condShort]); % SHORT labels
else
    conds = "Random";
end

rAll = vertcat(R.r);
rAll = rAll(isfinite(rAll));
rbins = linspace(min(rAll), max(rAll), RAD.nRBins+1);
rcent = (rbins(1:end-1) + rbins(2:end)) / 2;
nBins = numel(rcent);

for ic = 1:numel(conds)
    condS = conds(ic);

    if string(startType) == "csv"
        Rc = R([R.condShort] == condS);
    else
        Rc = R;
    end
    if isempty(Rc), continue; end

    medB = nan(nBins, numel(Rc)); loB = nan(nBins, numel(Rc)); hiB = nan(nBins, numel(Rc));
    medA = nan(nBins, numel(Rc)); loA = nan(nBins, numel(Rc)); hiA = nan(nBins, numel(Rc));

    for k = 1:numel(Rc)
        rFull = Rc(k).r(:);

        [eBfull, ~] = pickMetric(Rc(k).base.stats, METRIC);
        [eAfull, ~] = pickMetric(Rc(k).aft.stats,  METRIC);

        mB = true(size(rFull));
        mA = true(size(rFull));
        if isfield(Rc(k),'maskBase'), mB = Rc(k).maskBase(:); end
        if isfield(Rc(k),'maskAft'),  mA = Rc(k).maskAft(:);  end

        rB = rFull(mB); eB = eBfull(mB);
        rA = rFull(mA); eA = eAfull(mA);

        okB = isfinite(rB) & isfinite(eB);
        okA = isfinite(rA) & isfinite(eA);
        rB = rB(okB); eB = eB(okB);
        rA = rA(okA); eA = eA(okA);

        if isfinite(cap)
            keepB = eB <= cap;
            keepA = eA <= cap;
            rB = rB(keepB); eB = eB(keepB);
            rA = rA(keepA); eA = eA(keepA);
        end

        [mB_, loB_, hiB_] = radialQuantiles(rB, eB, rbins, RAD.qLo, RAD.qHi);
        [mA_, loA_, hiA_] = radialQuantiles(rA, eA, rbins, RAD.qLo, RAD.qHi);

        medB(:,k) = mB_; loB(:,k) = loB_; hiB(:,k) = hiB_;
        medA(:,k) = mA_; loA(:,k) = loA_; hiA(:,k) = hiA_;
    end

    mB = median(medB,2,'omitnan'); lB = median(loB,2,'omitnan'); hB = median(hiB,2,'omitnan');
    mA = median(medA,2,'omitnan'); lA = median(loA,2,'omitnan'); hA = median(hiA,2,'omitnan');

    if string(startType) == "csv" && N == 4 && ...
            (condS == "Tetrahedron" || contains(lower(condS), "tetra"))
        cBefore = HIST.cTetBefore;
        cAfter  = HIST.cTetAfter;
    else
        cBefore = HIST.cPyrBefore;
        cAfter  = HIST.cPyrAfter;
    end

    fillBand(ax, rcent, lB, hB, cBefore, HIST.alphaBefore);
    fillBand(ax, rcent, lA, hA, cAfter,  HIST.alphaAfter);

    plot(ax, rcent, mB, '--', 'Color', cBefore, 'LineWidth', RAD.beforeLineWidth);
    plot(ax, rcent, mA, '-',  'Color', cAfter,  'LineWidth', RAD.afterLineWidth);
end

if METRIC == 'pos'
    yline(ax, thr, '--k', 'thr', 'LineWidth', RAD.thresholdLineWidth);
end
xlabel(ax,'Distance from array centre (m)','Interpreter','latex');
ylabel(ax, pickMetricName(METRIC), 'Interpreter','latex');
xlim(ax, [rbins(1) rbins(end)]);

% if string(startType) == "csv"
%     L = strings(0);
%     for ic = 1:numel(conds)
%         L(end+1) = conds(ic) + " BEFORE"; %#ok<AGROW>
%         L(end+1) = conds(ic) + " AFTER";  %#ok<AGROW>
%     end
%     L(end+1) = "thr";
%     legend(ax, cellstr(L), 'Interpreter','latex', 'Location','best');
% else
%     legend(ax, ...
%         {sprintf('before %d--%d\\%% band', RAD.qLo, RAD.qHi), sprintf('after %d--%d\\%% band', RAD.qLo, RAD.qHi), ...
%          'before median','after median','thr'}, ...
%         'Interpreter','latex', 'Location','best');
% end

formatLatex(ax);
hold(ax,'off');
end

function HIST = ensureHistogramStyle(HIST)
% Supply shared defaults when a plotting section is run with an older
% HIST structure retained in the MATLAB workspace.
if ~isfield(HIST,'alphaBefore'), HIST.alphaBefore = 0.75; end
if ~isfield(HIST,'alphaAfter'),  HIST.alphaAfter  = 0.80; end
if ~isfield(HIST,'cPyrBefore'),  HIST.cPyrBefore  = [0.35 0.35 0.35]; end
if ~isfield(HIST,'cPyrAfter'),   HIST.cPyrAfter   = [0.05 0.38 0.78]; end
if ~isfield(HIST,'cTetAfter'),   HIST.cTetAfter   = [0.90 0.32 0.04]; end
if ~isfield(HIST,'cTetBefore')
    HIST.cTetBefore = (1-0.42)*HIST.cTetAfter + 0.42*[1 1 1];
end
end

function [med, qlo, qhi] = radialQuantiles(r, e, rbins, qLo, qHi)
r = r(:); e = e(:);
if numel(r) ~= numel(e)
    error('radialQuantiles:SizeMismatch', 'r (%d) and e (%d) must match.', numel(r), numel(e));
end

ok = isfinite(r) & isfinite(e);
r = r(ok); e = e(ok);

n = numel(rbins)-1;
med = nan(n,1); qlo = nan(n,1); qhi = nan(n,1);

for i = 1:n
    m = (r >= rbins(i)) & (r < rbins(i+1));
    if any(m)
        med(i) = median(e(m));
        qlo(i) = prctile(e(m), qLo);
        qhi(i) = prctile(e(m), qHi);
    end
end
end

function fillBand(ax, x, lo, hi, col, alpha)
x = x(:); lo = lo(:); hi = hi(:);
ok = isfinite(x) & isfinite(lo) & isfinite(hi);
x = x(ok); lo = lo(ok); hi = hi(ok);
if numel(x) < 2, return; end

fill(ax, [x; flipud(x)], [lo; flipud(hi)], col, ...
    'FaceAlpha', alpha, 'EdgeColor','none');
end

function run = pickRepresentativeRun(csvRuns, cond)
sel = arrayfun(@(r) r.condition==cond, csvRuns);
R = csvRuns(sel);

if isempty(R)
    error('pickRepresentativeRun:NoRuns','No runs found for condition "%s".', cond);
end

passA = arrayfun(@(r) r.passAfter, R);
p95A  = arrayfun(@(r) r.aft.stats.p95Err, R);

bestPass = max(passA);
cand = find(passA == bestPass);

if numel(cand) == 1
    run = R(cand);
else
    [~, j] = min(p95A(cand));
    run = R(cand(j));
end
end

function condsOut = orderConditions(condsIn, preferOrder)
condsIn = condsIn(:);
used = false(size(condsIn));

condsOut = strings(0,1);
for p = 1:numel(preferOrder)
    hit = contains(lower(condsIn), lower(preferOrder(p))) & ~used;
    if any(hit)
        condsOut = [condsOut; condsIn(hit)]; %#ok<AGROW>
        used(hit) = true;
    end
end

rest = condsIn(~used);
rest = sort(rest);
condsOut = [condsOut; rest];
end

function plot3D(ax, P, M, gridPts, vals, usedMask, labelStr, thr, cmax, catThr)
cla(ax); hold(ax,'on'); grid(ax,'on'); axis(ax,'equal'); view(ax,3);
xlabel(ax,'$x$ (m)','Interpreter','latex');
ylabel(ax,'$y$ (m)','Interpreter','latex');
zlabel(ax,'$z$ (m)','Interpreter','latex');

[sx,sy,sz] = sphere(30);
surf(ax, P.Lmax*sx, P.Lmax*sy, P.Lmax*sz, 'FaceAlpha',0.05,'EdgeAlpha',0.12);
surf(ax, P.R*sx,    P.R*sy,    P.R*sz,    'FaceAlpha',0.03,'EdgeAlpha',0.08);

plot3(ax,0,0,0,'pentagram','MarkerSize',6,'LineWidth',2,...
    'MarkerFaceColor','r','MarkerEdgeColor','none');

cmap = lines(size(M,1));
for ii = 1:size(M,1)
    scatter3(ax, M(ii,1),M(ii,2),M(ii,3), 40, cmap(ii,:), 'filled', ...
        'MarkerEdgeColor','k','LineWidth',1.5);
    text(ax, M(ii,1),M(ii,2),M(ii,3), sprintf('  %d',ii), ...
        'FontWeight','bold','FontSize',10,'Color',cmap(ii,:), 'Interpreter','none');
end

pts = gridPts(usedMask,:);
vv  = vals(usedMask);

exceedMask = vv > cmax;
catMask    = vv > catThr;

vvCap = min(vv, cmax);
scatter3(ax, pts(:,1),pts(:,2),pts(:,3), 28, vvCap, 'filled', 'MarkerFaceAlpha',0.75);

% if any(exceedMask)
%     scatter3(ax, pts(exceedMask,1),pts(exceedMask,2),pts(exceedMask,3), 20, 'm', 'filled', 'MarkerEdgeColor', 'none');
% end
% if any(catMask)
%     scatter3(ax, pts(catMask,1),pts(catMask,2),pts(catMask,3), 20, 'k', 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha',0.25);
% end

colormap(ax, turbo);
cb = colorbar(ax);
cb.Label.String = labelStr;
cb.Label.Interpreter = 'latex';
cb.TickLabelInterpreter = 'latex';
caxis(ax,[0, cmax]);

% text(ax, 0.02, 0.98, sprintf('thr=%.1f | 2\\times thr=%.1f | cat=%.1f', thr, cmax, catThr), ...
%     'Units','normalized', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
%     'FontWeight','bold', 'Interpreter','latex', ...
%     'BackgroundColor',[1 1 1 0.6], 'Margin',4);

formatLatex(ax);
hold(ax,'off');
end

function [vB, vA, thr] = pooledBeforeAfterVectors(R, METRIC)
thr = median(arrayfun(@(r) r.P.threshold_cm, R), 'omitnan');

vB = [];
vA = [];

for k = 1:numel(R)
    [eBfull, ~] = pickMetric(R(k).base.stats, METRIC);
    [eAfull, ~] = pickMetric(R(k).aft.stats,  METRIC);

    mB = true(size(eBfull));
    mA = true(size(eAfull));
    if isfield(R(k),'maskBase'), mB = R(k).maskBase(:); end
    if isfield(R(k),'maskAft'),  mA = R(k).maskAft(:);  end

    eB = eBfull(mB);
    eA = eAfull(mA);

    eB = eB(isfinite(eB));
    eA = eA(isfinite(eA));

    vB = [vB; eB(:)]; %#ok<AGROW>
    vA = [vA; eA(:)]; %#ok<AGROW>
end
end

function D = makeDeltaTableFromRuns(Runs, METRIC, useCap, capMult, COND_LABELS)
tmpl = struct( ...
    'startType', "", ...
    'condition', "", ...
    'N', NaN, ...
    'dPass_pp', NaN, ...
    'dMean', NaN, ...
    'dMedian', NaN, ...
    'dIQR', NaN, ...
    'dP95', NaN, ...
    'dTail_pp', NaN, ...
    'rowLabel', "" );

rows = repmat(tmpl, 0, 1);

csvRuns = Runs(arrayfun(@(x) x.startType=="csv", Runs));
if ~isempty(csvRuns)
    Ns = unique([csvRuns.N]);
    for iN = 1:numel(Ns)
        N = Ns(iN);
        Rn = csvRuns([csvRuns.N] == N);
        conds = unique([Rn.condition]);

        for ic = 1:numel(conds)
            cond = conds(ic);
            Rc = Rn([Rn.condition] == cond);
            if isempty(Rc), continue; end

            d = computeDeltasForGroup(Rc, METRIC, useCap, capMult, "csv", cond, N);
            condS = condLabel(cond, COND_LABELS);
            d.rowLabel = sprintf('%s', escapeLatex(condS));
            rows(end+1) = d; %#ok<AGROW>
        end
    end
end

rndRuns = Runs(arrayfun(@(x) x.startType=="random", Runs));
if ~isempty(rndRuns)
    Ns = unique([rndRuns.N]);
    for iN = 1:numel(Ns)
        N = Ns(iN);
        Rn = rndRuns([rndRuns.N] == N);
        if isempty(Rn), continue; end

        d = computeDeltasForGroup(Rn, METRIC, useCap, capMult, "random", "randomStart", N);
        d.rowLabel = sprintf('Random, $N_{mic}=%d$', N);
        rows(end+1) = d; %#ok<AGROW>
    end
end

D = struct2table(rows);
D.startType = string(D.startType);
D.condition = string(D.condition);
D.rowLabel  = string(D.rowLabel);
D = sortrows(D, {'startType','N','condition'});
end

function out = computeDeltasForGroup(Rgroup, METRIC, useCap, capMult, startType, condition, N)
[vB, vA, thr] = pooledBeforeAfterVectors(Rgroup, METRIC);
cap = capMult * thr;

tailB = 100 * mean(vB > cap, 'omitnan');
tailA = 100 * mean(vA > cap, 'omitnan');

passB = 100 * mean(vB <= thr, 'omitnan');
passA = 100 * mean(vA <= thr, 'omitnan');

if useCap
    vBstat = min(vB, cap);
    vAstat = min(vA, cap);
else
    vBstat = vB;
    vAstat = vA;
end

meanB = mean(vBstat, 'omitnan'); meanA = mean(vAstat, 'omitnan');
medB  = median(vBstat, 'omitnan'); medA  = median(vAstat, 'omitnan');

qB = prctile(vBstat, [25 75]); iqrB = qB(2) - qB(1);
qA = prctile(vAstat, [25 75]); iqrA = qA(2) - qA(1);

p95B = prctile(vBstat, 95);
p95A = prctile(vAstat, 95);

out = struct();
out.startType = string(startType);
out.condition = string(condition);
out.N = N;

out.dPass_pp = passA - passB;
out.dMean    = meanA - meanB;
out.dMedian  = medA  - medB;
out.dIQR     = iqrA  - iqrB;
out.dP95     = p95A  - p95B;
out.dTail_pp = tailA - tailB;

out.rowLabel = "";
end

function cmap = divergingBlueRedGreen(n)
if nargin < 1, n = 256; end
x = linspace(-1,1,n)';

cmap = zeros(n,3);

blue = [0.12 0.30 0.65];
red  = [0.80 0.15 0.15];
green= [0.10 0.55 0.20];

for i = 1:n
    if x(i) < 0
        t = abs(x(i));
        cmap(i,:) = (1-t)*blue + t*red;
    else
        t = x(i);
        cmap(i,:) = (1-t)*blue + t*green;
    end
end
end
