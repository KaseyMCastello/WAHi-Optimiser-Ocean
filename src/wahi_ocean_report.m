function fig = wahi_ocean_report(H, AFTER)
% wahi_ocean_report  Single-run report for a WAH-i Ocean optimisation history.
%
%   fig = wahi_ocean_report(H)            % baseline vs final
%   fig = wahi_ocean_report(H, "best")    % baseline vs best pass rate
%
% H is the run history exported by the GUI ("Export run") with fields
% H.meta (P, env, gridPts, gridInfo, ...) and H.iter (it, M, stats, J, ...).
% A second figure lists hydrophone start/end coordinates.

if nargin < 2, AFTER = 'final'; end
P = H.meta.P;
env = H.meta.env;
G = H.meta.gridPts;
thr = P.threshold_m;

passAll = arrayfun(@(x) x.stats.passRate, H.iter(:))';
medAll  = arrayfun(@(x) x.stats.medianErr, H.iter(:))';
itVals  = arrayfun(@(x) x.it, H.iter(:))';
accepted = arrayfun(@(x) logical(x.accepted), H.iter(:))';
isShake = arrayfun(@(x) ~isempty(strfind(char(x.note), 'shake')), H.iter(:))';

iBase = 1;
if strcmpi(char(AFTER), 'best')
    [~, ord] = sortrows([-passAll(:) medAll(:)]);
    iAfter = ord(1);
else
    iAfter = numel(H.iter);
end
base = H.iter(iBase); aft = H.iter(iAfter);
A = struct('ox', P.ox, 'oy', P.oy, 'mode', P.placeMode, 'refDepth', P.refDepth);
for f = {'yaw','pitch','roll'}
    if isfield(P, f{1}), A.(f{1}) = P.(f{1}); end
end
Hb = HydrophoneGeometry.toAbsolute(base.M, A, env);
Ha = HydrophoneGeometry.toAbsolute(aft.M, A, env);

makeHydTable(base.M, aft.M, AFTER);

fig = figure('Name', 'WAH-i Ocean report', 'Color', [1 1 1], 'Position', [40 40 1760 1040]);
t = tiledlayout(fig, 3, 4, 'Padding', 'compact', 'TileSpacing', 'compact');

cmax = 2*thr;
ax = nexttile(t, 1);
planMap(ax, env, G, base.stats, Hb, cmax);
title(ax, sprintf('BEFORE (it %.1f): pass %.1f%%, median %.0f m', base.it, 100*base.stats.passRate, base.stats.medianErr));

ax = nexttile(t, 2);
planMap(ax, env, G, aft.stats, Ha, cmax);
title(ax, sprintf('AFTER (%s, it %.1f): pass %.1f%%, median %.0f m', AFTER, aft.it, 100*aft.stats.passRate, aft.stats.medianErr));

ax = nexttile(t, 3);
hold(ax, 'on'); grid(ax, 'on');
plot(ax, itVals, 100*passAll, 'LineWidth', 1.5);
lg = {'pass rate'};
if any(accepted), scatter(ax, itVals(accepted), 100*passAll(accepted), 40, 'filled'); lg{end+1} = 'accepted'; end
if any(isShake),  scatter(ax, itVals(isShake), 100*passAll(isShake), 55, 'd'); lg{end+1} = 'shake/reset'; end
xlabel(ax, 'Iteration'); ylabel(ax, 'Pass rate (%)'); title(ax, 'Pass-rate progression');
legend(ax, lg, 'Location', 'best'); hold(ax, 'off');

ax = nexttile(t, 4);
hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
scatter(ax, base.M(:,1), base.M(:,2), 90, [0.6 0.6 0.6], 'filled', 'Marker', 'v');
cm = lines(size(aft.M,1));
scatter(ax, aft.M(:,1), aft.M(:,2), 90, cm, 'filled', 'Marker', 'v', 'MarkerEdgeColor', 'k');
for q = 1:size(aft.M,1)
    plot(ax, [base.M(q,1) aft.M(q,1)], [base.M(q,2) aft.M(q,2)], ':', 'Color', cm(q,:));
    text(ax, aft.M(q,1), aft.M(q,2), sprintf('  H%d', q), 'Color', cm(q,:)*0.8, 'FontWeight', 'bold');
end
xlabel(ax, 'x rel. (m)'); ylabel(ax, 'y rel. (m)'); title(ax, 'Geometry: before (grey) vs after');
hold(ax, 'off');

ub = base.stats.usedMask; ua = aft.stats.usedMask;
ax = nexttile(t, 5);
distPlot(ax, base.stats.score_m(ub), aft.stats.score_m(ua), thr, 'Score error (m)');
ax = nexttile(t, 6);
bb = base.stats.brg_deg(ub & base.stats.flag == 1); ba = aft.stats.brg_deg(ua & aft.stats.flag == 1);
thrB = max(HydrophoneGeometry.prctileLocal([bb; ba], 50), 1e-3);
distPlot(ax, bb, ba, thrB, 'Bearing error (deg)');

cB = mean(Hb(:,1:2), 1); cA = mean(Ha(:,1:2), 1);
rB = sqrt(sum((G(:,1:2) - cB).^2, 2)); rA = sqrt(sum((G(:,1:2) - cA).^2, 2));
ax = nexttile(t, 7);
vsRange(ax, rB(ub), base.stats.score_m(ub), rA(ua), aft.stats.score_m(ua), thr, 'Score error (m)');
ax = nexttile(t, 8);
okB = ub & base.stats.flag == 1; okA = ua & aft.stats.flag == 1;
vsRange(ax, rB(okB), base.stats.brg_deg(okB), rA(okA), aft.stats.brg_deg(okA), thrB, 'Bearing error (deg)');

nI = numel(H.iter); nH = size(base.M,1);
dx = nan(nI, nH); dy = nan(nI, nH);
for i = 1:nI
    if size(H.iter(i).M,1) == nH
        dx(i,:) = (H.iter(i).M(:,1) - base.M(:,1))';
        dy(i,:) = (H.iter(i).M(:,2) - base.M(:,2))';
    end
end
ax = nexttile(t, 9);  driftPlot(ax, itVals, dx, '\Delta x (m)', 'Hydrophone drift: x');
ax = nexttile(t, 10); driftPlot(ax, itVals, dy, '\Delta y (m)', 'Hydrophone drift: y');

ax = nexttile(t, 11);
hold(ax, 'on'); grid(ax, 'on');
scatter(ax, rA(ua)/1e3, aft.stats.minSNR(ua), 12, [0.1 0.5 0.9], 'filled', 'MarkerFaceAlpha', 0.4);
yline(ax, P.DT, '--', 'DT');
xlabel(ax, 'Range from array (km)'); ylabel(ax, 'Min. received SNR (dB)');
title(ax, 'Received SNR (after)'); hold(ax, 'off');

ax = nexttile(t, 12);
axis(ax, 'off');
txt = sprintf(['Run summary\n-----------\n' ...
    'Hydrophones: %d   sources: %d\n' ...
    'Field: R = %.1f km, spacing %.0f m, depths %s m\n' ...
    'Threshold: %.0f m (%s error)\n' ...
    'Pass (before/after): %.1f%% / %.1f%%\n' ...
    'Median (before/after): %.1f / %.1f m\n' ...
    'P95 (before/after): %.1f / %.1f m\n' ...
    'Detected (after): %.1f%%\n' ...
    'Solve: %s, TDOA: %s\n' ...
    'Truth SSP: %s, localiser: %s\n' ...
    'Median |move|: %.1f m'], ...
    nH, size(G,1), P.R/1e3, P.step, mat2str(P.depths), thr, P.scoreMetric, ...
    100*base.stats.passRate, 100*aft.stats.passRate, base.stats.medianErr, aft.stats.medianErr, ...
    base.stats.p95Err, aft.stats.p95Err, 100*aft.stats.detRate, P.solveMode, P.tdoaMethod, ...
    env.ssp.type, P.locSSPMode, median(sqrt(sum((aft.M(:,1:2) - base.M(:,1:2)).^2, 2))));
text(ax, 0, 1, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontName', 'Consolas', 'FontSize', 10);

fig.UserData = struct('base', base, 'after', aft, 'P', P);
end

%% ---------------------------------------------------------------------
function planMap(ax, env, G, st, Hn, cmax)
hold(ax, 'on'); grid(ax, 'on');
if ~isempty(env.bathy)
    B = env.bathy;
    sx = max(1, ceil(numel(B.x)/150)); sy = max(1, ceil(numel(B.y)/150));
    D = B.depth(1:sy:end, 1:sx:end);
    if numel(unique(round(D(:)))) > 1
        contour(ax, B.x(1:sx:end)/1e3, B.y(1:sy:end)/1e3, D, 8, 'LineColor', [0.75 0.75 0.78]);
    end
end
u = st.usedMask;
ok = u & isfinite(st.score_m);
scatter(ax, G(ok,1)/1e3, G(ok,2)/1e3, 18, min(st.score_m(ok), cmax), 'filled');
f = st.flag >= 2;
if any(f), plot(ax, G(f,1)/1e3, G(f,2)/1e3, 'x', 'Color', [0.5 0.5 0.5]); end
scatter(ax, Hn(:,1)/1e3, Hn(:,2)/1e3, 50, 'k', 'filled', 'Marker', 'v');
colormap(ax, turbo); caxis(ax, [0 cmax]);
cb = colorbar(ax); cb.Label.String = 'Score error (m)';
axis(ax, 'equal'); xlabel(ax, 'East (km)'); ylabel(ax, 'North (km)');
hold(ax, 'off');
end

function distPlot(ax, vB, vA, thr, lab)
hold(ax, 'on'); grid(ax, 'on');
vB = vB(isfinite(vB)); vA = vA(isfinite(vA));
cut = 4*thr;
edges = linspace(0, cut, 40);
histogram(ax, vB(vB <= cut), edges, 'Normalization', 'probability', 'FaceAlpha', 0.35);
histogram(ax, vA(vA <= cut), edges, 'Normalization', 'probability', 'FaceAlpha', 0.35);
xline(ax, thr, '--k', 'threshold');
xlabel(ax, lab); ylabel(ax, 'Probability');
title(ax, sprintf('%s (<= %.3g)', lab, cut));
legend(ax, {'before', 'after'}, 'Location', 'best');
txt = sprintf('> 4 x threshold\nbefore: %.1f%%\nafter: %.1f%%', 100*mean(vB > cut), 100*mean(vA > cut));
text(ax, 0.98, 0.7, txt, 'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontSize', 9, ...
    'BackgroundColor', [1 1 1], 'EdgeColor', [0.7 0.7 0.7]);
hold(ax, 'off');
end

function vsRange(ax, rB, eB, rA, eA, thr, lab)
hold(ax, 'on'); grid(ax, 'on');
cap = 3*thr;
h1 = scatter(ax, rB/1e3, min(eB, cap), 10, [0.6 0.6 0.6], 'filled', 'MarkerFaceAlpha', 0.3);
h2 = scatter(ax, rA/1e3, min(eA, cap), 10, [0.1 0.5 0.9], 'filled', 'MarkerFaceAlpha', 0.4);
hh = [h1 h2]; lab = {'before', 'after'};
rr = [rB(:); rA(:)];
if numel(rr) > 2 && max(rr) > min(rr)
    edges = linspace(min(rr), max(rr), 12);
    mB = binMedian(rB, eB, edges); mA = binMedian(rA, eA, edges);
    ctr = (edges(1:end-1) + edges(2:end))/2/1e3;
    h3 = plot(ax, ctr, min(mB, cap), 'k--', 'LineWidth', 1.5);
    h4 = plot(ax, ctr, min(mA, cap), 'b-', 'LineWidth', 2);
    hh = [hh h3 h4]; lab = [lab {'before median', 'after median'}];
end
yline(ax, thr, '--', 'threshold');
xlabel(ax, 'Range from array centroid (km)'); ylabel(ax, lab);
title(ax, [lab ' vs range']);
legend(ax, hh, lab, 'Location', 'northwest');
hold(ax, 'off');
end

function m = binMedian(r, e, edges)
m = nan(1, numel(edges)-1);
for i = 1:numel(edges)-1
    s = r >= edges(i) & r <= edges(i+1) & isfinite(e);
    if any(s), m(i) = median(e(s)); end
end
end

function driftPlot(ax, it, D, ylab, ttl)
hold(ax, 'on'); grid(ax, 'on');
for m = 1:size(D,2), plot(ax, it, D(:,m), 'LineWidth', 1.2); end
xlabel(ax, 'Iteration'); ylabel(ax, ylab); title(ax, ttl);
legend(ax, compose('H%d', 1:size(D,2)), 'Location', 'best', 'Box', 'off');
hold(ax, 'off');
end

function makeHydTable(M0, M1, AFTER)
if size(M0,1) ~= size(M1,1), return; end
n = size(M0,1);
d = M1 - M0;
data = [(1:n)' M0 M1 d sqrt(sum(d.^2,2))];
f = figure('Name', sprintf('Hydrophone start/end (%s)', AFTER), 'Color', [1 1 1], 'Position', [120 120 1000 320]);
uitable(f, 'Data', data, 'ColumnName', {'H','x0','y0','z0','x1','y1','z1','dx','dy','dz','|d|'}, ...
    'Units', 'normalized', 'Position', [0.01 0.01 0.98 0.98]);
end
