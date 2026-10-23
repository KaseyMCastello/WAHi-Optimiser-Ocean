function [stats, diag] = wahi_evaluate_field(loc, gridPts, E)
% wahi_evaluate_field  Run simulate -> TDOA -> solve at every grid point.
%
%   [stats, diag] = wahi_evaluate_field(loc, gridPts, E)
%
% loc     : WhaleCallLocaliser
% gridPts : K x 3 source positions (world frame)
% E fields
%   threshold_m      pass threshold on the score error (m)
%   failPenalty_m    error assigned to failed / undetected points (m)
%   minDistToAnyHyd  sources closer than this to a hydrophone are skipped (m)
%   scoreMetric      'horizontal' | '3d'
%   undetected       'fail' (count as failures) | 'exclude'
%   seed             base seed; point i uses rng(seed + i) so every geometry
%                    sees the same noise realisation (common random numbers)
%   useParallel      true to use parfor (Parallel Computing Toolbox)
%   progressFcn      optional @(i,n) -> true to abort (serial only)
%
% Flag codes in stats.flag: 0 skipped (too close), 1 ok, 2 undetected,
%                           3 TDOA failure, 4 solver failure

if ~isfield(E,'scoreMetric'), E.scoreMetric = 'horizontal'; end
if ~isfield(E,'undetected'), E.undetected = 'fail'; end
if ~isfield(E,'seed') || ~isfinite(E.seed), E.seed = 1; end
if ~isfield(E,'useParallel'), E.useParallel = false; end
if ~isfield(E,'progressFcn'), E.progressFcn = []; end

nP = size(gridPts,1);
errH = nan(nP,1); err3 = nan(nP,1); errZ = nan(nP,1); errB = nan(nP,1); errR = nan(nP,1);
nDet = zeros(nP,1); minSNR = nan(nP,1); maxSNR = nan(nP,1);
flag = zeros(nP,1);
est = nan(nP,3);

Hn = loc.Hnom;
dmin = zeros(nP,1);
for i = 1:nP
    dmin(i) = min(sqrt(sum((Hn - gridPts(i,:)).^2,2)));
end
tooCloseMask = dmin < E.minDistToAnyHyd;

aborted = false; abortAt = nP;
seed0 = E.seed;
if E.useParallel && canUseParallel()
    parfor i = 1:nP
        if tooCloseMask(i), continue; end
        rng(seed0 + i, 'twister');
        o = loc.evaluatePoint(gridPts(i,:));
        [flag(i), errH(i), err3(i), errZ(i), errB(i), errR(i), nDet(i), minSNR(i), maxSNR(i)] = unpack(o);
        est(i,:) = o.est;
    end
else
    for i = 1:nP
        if tooCloseMask(i), continue; end
        rng(seed0 + i, 'twister');
        try
            o = loc.evaluatePoint(gridPts(i,:));
        catch
            o = WhaleCallLocaliser.emptyResult(); o.flag = 'solve';
        end
        [flag(i), errH(i), err3(i), errZ(i), errB(i), errR(i), nDet(i), minSNR(i), maxSNR(i)] = unpack(o);
        est(i,:) = o.est;
        if ~isempty(E.progressFcn) && mod(i, 10) == 0
            if E.progressFcn(i, nP), aborted = true; abortAt = i; break; end
        end
    end
end

usedMask = ~tooCloseMask;
if aborted, usedMask(abortAt+1:end) = false; end
if strcmpi(E.undetected, 'exclude')
    usedMask = usedMask & flag ~= 2;
end

if strcmpi(E.scoreMetric, '3d'), score = err3; else, score = errH; end
score(usedMask & ~isfinite(score)) = E.failPenalty_m;
score = min(score, E.failPenalty_m);
score(~usedMask) = NaN;

passMask = usedMask & score <= E.threshold_m;
nUsed = max(1, sum(usedMask));
u = score(usedMask);
if isempty(u), u = E.failPenalty_m; end

stats = struct();
stats.score_m = score;
stats.errH_m = errH; stats.err3_m = err3; stats.errZ_m = errZ;
stats.brg_deg = errB; stats.range_m = errR;
stats.est = est;
stats.nDet = nDet; stats.minSNR = minSNR; stats.maxSNR = maxSNR;
stats.flag = flag;
stats.usedMask = usedMask;
stats.passMask = passMask;
stats.passRate = sum(passMask)/nUsed;
stats.meanErr = mean(u);
stats.medianErr = median(u);
stats.p95Err = HydrophoneGeometry.prctileLocal(u, 95);
stats.maxErr = max(u);
stats.threshold_m = E.threshold_m;
stats.scoreMetric = E.scoreMetric;
stats.detRate = sum(flag == 1 | flag >= 3) / max(1, sum(~tooCloseMask));

diag = struct();
diag.nP = nP;
diag.gridUsed = sum(usedMask);
diag.tooClose = sum(tooCloseMask);
diag.undetected = sum(flag == 2);
diag.tdoaFail = sum(flag == 3);
diag.solveFail = sum(flag == 4);
diag.aborted = aborted;
end

function [f, eH, e3, eZ, eB, eR, nd, s1, s2] = unpack(o)
switch o.flag
    case 'ok',         f = 1;
    case 'undetected', f = 2;
    case 'tdoa',       f = 3;
    otherwise,         f = 4;
end
eH = o.errH; e3 = o.err3; eZ = o.errZ; eB = o.errBrg; eR = o.errRange;
nd = o.nDet; s1 = o.minSNR; s2 = o.maxSNR;
end

function ok = canUseParallel()
ok = false;
try
    ok = license('test','Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
    if ok && isempty(gcp('nocreate')), parpool; end
catch
    ok = false;
end
end
