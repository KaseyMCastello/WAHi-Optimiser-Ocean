function [H, runInfo] = iwah_run_optimiser_core(P, env, call, M0, seed, doOptimise)
% iwah_run_optimiser_core  Headless WAH-i Ocean evaluation / optimisation.
%
%   [H, runInfo] = iwah_run_optimiser_core(P, env, call, M0)
%   [H, runInfo] = iwah_run_optimiser_core(P, env, call, M0, seed, doOptimise)
%
% Runs the same pipeline as the GUI without any plotting, for batch studies.
% The easiest way to get the inputs is from the GUI: click "Export workspace"
% (MATLAB only), which creates wahi_P, wahi_env, wahi_call and wahi_M.
%
%   P          parameter struct as produced by the GUI (wahi_P)
%   env        OceanEnvironment struct incl. the truth sound-speed model (wahi_env)
%   call       call waveform at P.fs (wahi_call)
%   M0         N x 3 relative hydrophone geometry, or the path of a geometry CSV
%   seed       optional seed (default P.seed)
%   doOptimise true (default) to optimise, false to evaluate M0 only
%
% H has the same layout as the GUI's "Export run" history (H.meta, H.iter),
% so wahi_ocean_report(H) works on it.
%
% Example
%   load('my_session_inputs.mat')     % wahi_P, wahi_env, wahi_call
%   [H, info] = iwah_run_optimiser_core(wahi_P, wahi_env, wahi_call, 'configs/hydrophones_4ch_triangle_centre.csv');
%   wahi_ocean_report(H);

if nargin < 5 || isempty(seed), seed = P.seed; end
if nargin < 6 || isempty(doOptimise), doOptimise = true; end
if ischar(M0) || isstring(M0)
    csvPath = char(M0);
    M0 = HydrophoneGeometry.readCSV(csvPath);
else
    csvPath = '';
end
P.seed = seed;
rng(seed, 'twister');
env = OceanEnvironment.rebuild(env);

[fLo, fHi] = WhaleCallLibrary.estimateBand(call, P.fs);
F = struct('cx', P.fcx, 'cy', P.fcy, 'R', P.R, 'Rin', P.Rin, 'step', P.step, ...
    'layout', P.layout, 'azFrom', P.azFrom, 'azTo', P.azTo, 'depths', P.depths);
[grid, gridInfo] = wahi_build_field_grid(F, env);
if isempty(grid), error('No synthetic sources in the field.'); end

cfg = struct('fs', P.fs, 'call', call(:), 'band', [fLo fHi], 'snrMode', P.snrMode, ...
    'snrFixed', P.snrFixed, 'SL', P.SL, 'NL', P.NL, 'DT', P.DT, 'minDetect', P.minDetect, ...
    'timingSigma', P.timingSigma, 'posSigma', P.posSigma, 'solveMode', P.solveMode, ...
    'assumedDepth', P.assumedDepth, 'tdoaMethod', P.tdoaMethod, 'seed', seed, ...
    'locSSP', OceanEnvironment.constantSSP(P.cLoc));
if isfield(P, 'useTrueDepth'), cfg.useTrueDepth = P.useTrueDepth; end
E = struct('threshold_m', P.threshold_m, 'failPenalty_m', P.failPen, 'minDistToAnyHyd', P.minDistAny, ...
    'scoreMetric', P.scoreMetric, 'undetected', P.undetected, 'seed', seed, ...
    'useParallel', isfield(P,'useParallel') && P.useParallel);

evalM = @(M) evaluateGeometry(M, P, env, cfg, grid, E);
J = @(st) st.passRate - P.lambda*min(st.meanErr/P.failPen, 1);

M = M0;
[bestStats, bestDiag] = evalM(M);
bestJ = J(bestStats);

H = struct();
H.meta = struct('version', 'WAH-i Ocean 1.0', 'rngSeed', seed, 'timeStart', datestr(now), ...
    'P', P, 'env', env, 'call', call(:), 'band', [fLo fHi], 'gridPts', grid, 'gridInfo', gridInfo);
H.meta.start = struct('source', 'batch', 'csvPath', csvPath, 'M', M0);
H.iter = struct('it', {}, 'M', {}, 'hydIdx', {}, 'stepSize', {}, 'dirsTried', {}, ...
    'stats', {}, 'diag', {}, 'J', {}, 'accepted', {}, 'note', {});
H.iter(end+1) = rec(0, M, NaN, NaN, [], bestStats, bestDiag, bestJ, true, 'baseline');

if doOptimise
    ang = [0 180 90 270 45 225 135 315 22.5 202.5 112.5 292.5 67.5 247.5 157.5 337.5];
    dirsAll = [cosd(ang)' sind(ang)' zeros(numel(ang),1)];
    if P.moveZ, dirsAll = [dirsAll(1:4,:); 0 0 1; 0 0 -1; dirsAll(5:end,:)]; end
    N = size(M,1);
    stall = 0;
    for it = 1:P.maxIters
        if P.maxIters > 1
            stepFrac = max(P.minStep, P.initStep*(1 - (it-1)/(P.maxIters-1)));
        else
            stepFrac = P.initStep;
        end
        stepSize = stepFrac*P.Lmax;
        h = mod(it-1, N) + 1;
        dirs = dirsAll(1:min(P.probes, size(dirsAll,1)), :);
        rv = [M(h,1:2) - mean(M(:,1:2),1) 0];
        if norm(rv) > 1e-9, rv = rv/norm(rv); else, rv = [1 0 0]; end
        dirs = [dirs; rv; -rv]; %#ok<AGROW>

        bJ = -inf; bM = []; bS = []; bD = [];
        for k = 1:size(dirs,1)
            Mt = M;
            Mt(h,:) = Mt(h,:) + stepSize*dirs(k,:);
            Mt = HydrophoneGeometry.clampXY(HydrophoneGeometry.recentreXY(Mt), P.Lmax);
            if strcmp(P.placeMode, 'seafloor'), Mt(:,3) = max(Mt(:,3), 0); end
            if ~HydrophoneGeometry.isValid(Mt, P.Lmax, P.minPair), continue; end
            try
                [st, dg] = evalM(Mt);
            catch
                continue;
            end
            if J(st) > bJ, bJ = J(st); bM = Mt; bS = st; bD = dg; end
        end
        accepted = false;
        if ~isempty(bM) && bJ > bestJ + P.epsJ
            M = bM; bestJ = bJ; bestStats = bS; bestDiag = bD; stall = 0; accepted = true;
            note = 'accepted (improved J)';
        else
            stall = stall + 1; note = sprintf('rejected (stall=%d)', stall);
        end
        H.iter(end+1) = rec(it, M, h, stepSize, dirs, bestStats, bestDiag, bestJ, accepted, note); %#ok<AGROW>
        if bestStats.passRate >= P.targetPass, break; end
        if stall >= P.stall
            Ms = HydrophoneGeometry.shake(M, P.shake*P.Lmax, P.Lmax, P.minPair, P.moveZ);
            try
                [st, dg] = evalM(Ms);
                M = Ms; bestStats = st; bestDiag = dg; bestJ = J(st);
                H.iter(end+1) = rec(it + 0.5, M, NaN, NaN, [], bestStats, bestDiag, bestJ, true, 'shake/reset'); %#ok<AGROW>
            catch
            end
            stall = 0;
        end
    end
end

H.meta.timeEnd = datestr(now);
H.meta.final = struct('M', M, 'stats', bestStats, 'diag', bestDiag, 'J', bestJ);
runInfo = struct('seed', seed, 'N', size(M,1), 'csvPath', csvPath, 'mode', 'batch', ...
    'finalPass', bestStats.passRate, 'finalMedian', bestStats.medianErr);
end

function r = rec(it, M, h, stepSize, dirs, st, dg, J, accepted, note)
r = struct('it', it, 'M', M, 'hydIdx', h, 'stepSize', stepSize, 'dirsTried', dirs, ...
    'stats', st, 'diag', dg, 'J', J, 'accepted', accepted, 'note', char(note));
end

function [st, dg] = evaluateGeometry(M, P, env, cfg, grid, E)
A = struct('ox', P.ox, 'oy', P.oy, 'mode', P.placeMode, 'refDepth', P.refDepth);
for f = {'yaw','pitch','roll'}
    if isfield(P, f{1}), A.(f{1}) = P.(f{1}); end
end
Hn = HydrophoneGeometry.toAbsolute(M, A, env);
if any(OceanEnvironment.waterDepth(env, Hn(:,1), Hn(:,2)) <= 0)
    error('A hydrophone is on land.');
end
if any(HydrophoneGeometry.belowSeabed(Hn, env)), error('A hydrophone is below the seabed.'); end
Xp = Hn(:,1:2) - mean(Hn(:,1:2), 1);
if rank(Xp, 1e-6*max(1, max(abs(Xp(:))))) < 2, error('Hydrophones are collinear in plan view.'); end
cfg = wahi_geometry_cfg(cfg, P, env, Hn, grid);
loc = WhaleCallLocaliser(cfg, env, Hn);
[st, dg] = wahi_evaluate_field(loc, grid, E);
end
