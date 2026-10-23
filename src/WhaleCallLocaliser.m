% WhaleCallLocaliser - simulate and localise a marine-mammal call on a
% hydrophone array (TDOA multilateration). Replaces BatCallLocaliser.
%
%   loc = WhaleCallLocaliser(cfg, env, Hnom)
%   out = loc.evaluatePoint(src)            % one source, full pipeline
%   sim = loc.simulate(src)                 % time-domain signals (for plots)
%
% Inputs
%   env  : OceanEnvironment struct (truth propagation: bathymetry + SSP)
%   Hnom : N x 3 nominal (surveyed) hydrophone positions, world frame
%   cfg  : struct with fields
%     fs, call (column vector), band ([fLo fHi] Hz)
%     snrMode ('sonar'|'fixed'), snrFixed, SL, NL (dB), DT (dB), minDetect
%     timingSigma (s), posSigma (m)  - clock and position uncertainty
%     solveMode ('2d'|'3d'), assumedDepth (m, 2-D solve)
%     useTrueDepth : 2-D solve uses each source's true depth (known depth)
%     tdoaMethod ('waveform'|'envelope')
%     locSSP  : sound speed model assumed by the localiser (constant/profile)
%     searchRmax (m) : outer radius for the initial coarse search
%     seed    : seed for the position-uncertainty realisation (optional)
%     propModel : 'straight' (default) or 'ray' (OceanRayModel: refraction,
%                 ducting, shadow zones, surface/bottom multipath)
%     ray     : struct for propModel 'ray': ssp (profile), D (flat bottom, m),
%               o (OceanRayModel options incl. BL), arrivals ('all' |
%               'strongest' | 'first'), window (s, multipath delay window)
%     locModel  : 'straight' (default, uses locSSP) or 'ray' (ray-traced
%                 first-arrival tables built from locRay.ssp / locRay.D)
%     locRay  : struct(ssp, D, o, nZ) for locModel 'ray'
%
% Signal model (per source)
%   - straight-ray travel times through env (incl. depth/3-D sound speed)
%   - received SNR per hydrophone from the passive sonar equation
%     SNR = SL - TL(r, f) - NL, or a fixed SNR
%   - independent band-limited Gaussian noise on every hydrophone
%   - only hydrophones with SNR >= DT take part; at least minDetect needed
%   - signals are built and cross-correlated in the frequency domain
%     (exact fractional delays), TDOA from the correlation peak with
%     parabolic sub-sample interpolation
%   - position solved by a coarse polar search + Levenberg-Marquardt using
%     the localiser's own sound-speed assumption and NOMINAL hydrophone
%     positions (so SSP mismatch, survey error and clock error all show up
%     as localisation error)

classdef WhaleCallLocaliser
    properties
        cfg
        env
        Hnom
        Htrue
        Nfft
        Xcall
        fvec
        bandMask
        Kb
        Ps
        maxLag
        maxDepth
        cRef
        fc
    end

    methods
        function obj = WhaleCallLocaliser(cfg, env, Hnom)
            obj.cfg = WhaleCallLocaliser.fillDefaults(cfg);
            obj.env = env;
            obj.Hnom = Hnom;
            cfg = obj.cfg;

            if isfield(cfg,'seed') && ~isempty(cfg.seed) && isfinite(cfg.seed)
                rng(cfg.seed, 'twister');
            end
            obj.Htrue = Hnom;
            if cfg.posSigma > 0
                obj.Htrue = Hnom + cfg.posSigma*randn(size(Hnom));
            end

            fs = cfg.fs;
            call = cfg.call(:);
            Lc = numel(call);
            D = HydrophoneGeometry.pairDist(Hnom);
            cmin = 1400;
            obj.maxLag = ceil((max(D(:))/cmin + 6*cfg.timingSigma)*fs) + 8;
            extra = 0;
            if strcmpi(cfg.propModel, 'ray') && strcmpi(cfg.ray.arrivals, 'all')
                extra = ceil(cfg.ray.window*fs);     % room for the multipath delay spread
            end
            obj.Nfft = 2^nextpow2(Lc + 2*obj.maxLag + extra + 32);
            obj.Xcall = fft(call, obj.Nfft);
            k = (0:obj.Nfft-1)';
            f = k*fs/obj.Nfft;
            f(k > obj.Nfft/2) = f(k > obj.Nfft/2) - fs;
            obj.fvec = f;
            band = cfg.band;
            obj.bandMask = abs(f) >= band(1) & abs(f) <= band(2);
            if ~any(obj.bandMask), obj.bandMask = true(size(f)); end
            obj.Kb = sum(obj.bandMask);
            obj.Ps = sum(call.^2)/max(Lc,1);
            obj.maxDepth = OceanEnvironment.maxWaterDepth(env);
            obj.cRef = cfg.locSSP.c0;
            obj.fc = mean(band);
        end

        %% ------------------------------------------------------------------
        function [snr, r, arr] = receivedSNR(obj, src, H)
            % Received SNR per hydrophone (1 x N). For the ray model also
            % returns the eigenrays (cell, one struct per hydrophone).
            r = sqrt(sum((H - src).^2, 2))';
            arr = {};
            if strcmpi(obj.cfg.propModel, 'ray')
                arr = obj.rayArrivals(src, H);
                TL = inf(size(r));
                for j = 1:numel(arr)
                    if ~isempty(arr{j}.t), TL(j) = -10*log10(sum(10.^(-arr{j}.TL/10))); end
                end
            else
                TL = OceanEnvironment.transmissionLoss(obj.env, r, obj.fc);
            end
            if strcmpi(obj.cfg.snrMode, 'fixed')
                snr = obj.cfg.snrFixed * ones(size(r));
                snr(~isfinite(TL)) = -inf;          % ray shadow zone: nothing arrives
            else
                snr = obj.cfg.SL - TL - obj.cfg.NL;
            end
        end

        function arr = rayArrivals(obj, src, H)
            R = obj.cfg.ray;
            arr = cell(1, size(H,1));
            zs = -src(3);
            for j = 1:size(H,1)
                S = WhaleCallLocaliser.cachedPaths(R.ssp, R.D, -H(j,3), zs, R.o);
                rh = norm(src(1:2) - H(j,1:2));
                arr{j} = OceanRayModel.eigenrays(S, rh, obj.fc, obj.env.prop.absorption);
            end
        end

        function [Y, t, snr, idx] = buildSpectra(obj, src)
            % Frequency-domain received signals for detecting hydrophones.
            % t = (first) arrival time per hydrophone, snr = SNR per hydrophone.
            isRay = strcmpi(obj.cfg.propModel, 'ray');
            N = size(obj.Htrue,1);
            if isRay
                [snr, ~, arr] = obj.receivedSNR(src, obj.Htrue);
                t = nan(1, N);
                for j = 1:N
                    A = arr{j};
                    if isempty(A.t), continue; end
                    switch lower(obj.cfg.ray.arrivals)
                        case 'first',     k = 1;
                        case 'strongest', [~, k] = min(A.TL);
                        otherwise,        k = find(A.TL <= min(A.TL) + 30);
                    end
                    if numel(k) < numel(A.t)      % only these paths exist in the signal
                        fn = {'t','TL','sgn','nb','ns','th0'};
                        for f = 1:numel(fn), A.(fn{f}) = A.(fn{f})(k); end
                        if ~strcmpi(obj.cfg.ray.arrivals, 'all') && ~strcmpi(obj.cfg.snrMode, 'fixed')
                            snr(j) = obj.cfg.SL - (-10*log10(sum(10.^(-A.TL/10)))) - obj.cfg.NL;
                        end
                        arr{j} = A;
                    end
                    t(j) = min(A.t);
                end
            else
                t = OceanEnvironment.travelTime(obj.env, src, obj.Htrue);    % 1 x N
                snr = obj.receivedSNR(src, obj.Htrue);
            end
            dtClk = zeros(1, N);
            if obj.cfg.timingSigma > 0
                dtClk = obj.cfg.timingSigma*randn(1, N);
            end
            t = t + dtClk;
            idx = find(snr >= obj.cfg.DT & isfinite(t));
            Y = [];
            if numel(idx) < obj.cfg.minDetect, return; end
            n = numel(idx);
            % arrival sets (times, amplitudes) per detecting hydrophone
            tt = cell(1, n); aa = cell(1, n);
            for q = 1:n
                j = idx(q);
                if isRay
                    A = arr{j};
                    tt{q} = A.t' + dtClk(j);
                    a = 10.^(-(A.TL' - min(A.TL))/20) .* A.sgn';
                    aa{q} = a / sqrt(sum(a.^2));
                else
                    tt{q} = t(j); aa{q} = 1;
                end
            end
            t0 = min(cellfun(@min, tt));
            maxShift = (obj.Nfft - numel(obj.cfg.call) - obj.maxLag) / obj.cfg.fs;
            Y = zeros(obj.Nfft, n);
            sc = sqrt(obj.Nfft / obj.Kb);
            for q = 1:n
                rel = tt{q} - t0;
                a = aa{q};
                keep = rel <= maxShift;
                if ~any(keep), keep = 1; rel = min(rel, maxShift); end
                H = exp(-1i*2*pi*obj.fvec*rel(keep)) * a(keep)';
                S = obj.Xcall .* H;
                sig2 = obj.Ps / 10^(snr(idx(q))/10);
                W = fft(randn(obj.Nfft,1)) * (sqrt(sig2)*sc);
                Y(:,q) = (S + W) .* obj.bandMask;
            end
        end

        function tdoa = estimateTDOA(obj, Y, ref)
            n = size(Y,2);
            useAbs = true;
            if strcmpi(obj.cfg.tdoaMethod, 'envelope')
                Y = WhaleCallLocaliser.envelopeSpectra(Y);
                useAbs = false;
            end
            L = obj.maxLag; N = obj.Nfft;
            tdoa = zeros(1,n);
            Yr = conj(Y(:,ref));
            for j = 1:n
                if j == ref, continue; end
                r = real(ifft(Y(:,j) .* Yr));
                w = [r(N-L+1:N); r(1:L+1)];      % lags -L..L
                if useAbs
                    [~, m] = max(abs(w));        % polarity may flip (surface reflection)
                    w = w * sign(w(m) + (w(m) == 0));
                else
                    [~, m] = max(w);
                end
                d = 0;
                if m > 1 && m < numel(w)
                    y1 = w(m-1); y2 = w(m); y3 = w(m+1);
                    den = y1 - 2*y2 + y3;
                    if den < 0, d = 0.5*(y1 - y3)/den; end
                end
                tdoa(j) = ((m - 1 - L) + d) / obj.cfg.fs;
            end
        end

        function out = evaluatePoint(obj, src)
            out = WhaleCallLocaliser.emptyResult();
            [Y, ~, snr, idx] = obj.buildSpectra(src);
            out.minSNR = max(min(snr), -99);
            out.maxSNR = max(max(snr), -99);
            out.nDet = numel(idx);
            if isempty(Y)
                out.flag = 'undetected';
                return;
            end
            [~, ref] = max(snr(idx));
            tdoa = obj.estimateTDOA(Y, ref);
            if any(~isfinite(tdoa)), out.flag = 'tdoa'; return; end
            est = obj.solve(tdoa, obj.Hnom(idx,:), ref, src(3));
            if any(~isfinite(est)), out.flag = 'solve'; return; end
            out = WhaleCallLocaliser.fillErrors(out, est, src, mean(obj.Hnom,1));
            out.flag = 'ok';
        end

        function sim = simulate(obj, src)
            % Time-domain signals and diagnostics for a single source (plots)
            [Y, t, snr, idx] = obj.buildSpectra(src);
            sim = struct('t',t,'snr',snr,'idx',idx,'fs',obj.cfg.fs,'signals',[],'tdoa',[],'est',NaN(1,3));
            if isempty(Y), return; end
            sim.signals = real(ifft(Y));
            [~, ref] = max(snr(idx));
            sim.ref = ref;
            sim.tdoa = obj.estimateTDOA(Y, ref);
            sim.trueTdoa = t(idx) - t(idx(ref));
            sim.est = obj.solve(sim.tdoa, obj.Hnom(idx,:), ref, src(3));
        end

        %% ------------------------------------------------------------------
        function est = solve(obj, tdoa, Hs, ref, zTrue)
            % zTrue (optional): true source elevation, used by the 2-D solver
            % when cfg.useTrueDepth is set (perfect depth knowledge)
            cfg = obj.cfg;
            obs = tdoa(:)';
            ctr = mean(Hs(:,1:2),1);
            is3 = strcmpi(cfg.solveMode, '3d');
            zFix = -abs(cfg.assumedDepth);
            if nargin >= 5 && cfg.useTrueDepth && isfinite(zTrue), zFix = zTrue; end
            zmin = -obj.maxDepth;

            % coarse polar search
            az = (0:10:350)*pi/180;
            rg = logspace(log10(10), log10(max(cfg.searchRmax, 100)), 26);
            [AA, RR] = meshgrid(az, rg);
            cx = ctr(1) + RR(:).*sin(AA(:));
            cy = ctr(2) + RR(:).*cos(AA(:));
            if is3
                dz = linspace(2, obj.maxDepth - 2, 6);
                nC = numel(cx);
                cx = repmat(cx, numel(dz), 1); cy = repmat(cy, numel(dz), 1);
                cz = -kron(dz(:), ones(nC,1));
            else
                cz = zFix*ones(size(cx));
            end
            C = [cx cy cz];
            TB = obj.locTables(Hs, zFix, is3);
            T = obj.locTimes(C, Hs, TB);
            Dt = T - T(:,ref);
            cost = sum((Dt - obs).^2, 2);
            [~, ord] = sort(cost);
            starts = C(ord(1),:);
            for k = 2:min(numel(ord), 400)
                c = C(ord(k),:);
                rr = max(norm(starts(1,1:2) - ctr), 100);
                if all(sqrt(sum((starts(:,1:2) - c(1:2)).^2, 2)) > 0.15*rr)
                    starts = [starts; c]; %#ok<AGROW>
                end
                if size(starts,1) >= 3, break; end
            end

            bestCost = inf; est = NaN(1,3);
            for s = 1:size(starts,1)
                if is3, u0 = starts(s,:); else, u0 = starts(s,1:2); end
                [u, cst] = obj.levmar(u0, obs, Hs, ref, is3, zFix, zmin, TB);
                if cst < bestCost
                    bestCost = cst;
                    if is3, est = u; else, est = [u zFix]; end
                end
            end
        end

        function [u, cst] = levmar(obj, u, obs, Hs, ref, is3, zFix, zmin, TB)
            c0 = obj.cRef;
            resid = @(uu) obj.residualM(uu, obs, Hs, ref, c0, is3, zFix, TB);
            r = resid(u); cst = r'*r;
            mu = 1e-3; d = numel(u);
            for it = 1:60
                h = max(0.5, 1e-6*norm(u));
                J = zeros(numel(r), d);
                for q = 1:d
                    up = u; up(q) = up(q) + h;
                    J(:,q) = (resid(up) - r)/h;
                end
                A = J'*J; g = J'*r;
                improved = false;
                for tries = 1:8
                    step = -(A + mu*diag(diag(A) + 1e-9)) \ g;
                    un = u + step';
                    if is3, un(3) = min(max(un(3), zmin), 0); end
                    rn = resid(un); cn = rn'*rn;
                    if cn < cst
                        u = un; r = rn; cst = cn; mu = max(mu/3, 1e-9); improved = true; break;
                    end
                    mu = mu*4;
                end
                if ~improved || norm(step) < 1e-3, break; end
            end
        end

        function r = residualM(obj, u, obs, Hs, ref, c0, is3, zFix, TB)
            if is3, C = u(:)'; else, C = [u(:)' zFix]; end
            T = obj.locTimes(C, Hs, TB);
            r = c0*((T - T(ref)) - obs)';
        end

        function TB = locTables(obj, Hs, zFix, is3)
            % Ray-traced first-arrival tables for the localiser (empty for straight rays)
            TB = [];
            if ~strcmpi(obj.cfg.locModel, 'ray'), return; end
            L = obj.cfg.locRay;
            rmax = 1.3*obj.cfg.searchRmax + 2000;
            if is3
                zl = linspace(1, L.D - 1, L.nZ);
            else
                zl = min(max(-zFix, 0.5), L.D - 0.5);
            end
            TB = struct('zl', zl, 'T', {cell(size(Hs,1), numel(zl))}, 'dr', NaN, 'rmax', rmax);
            for i = 1:size(Hs,1)
                for l = 1:numel(zl)
                    tb = WhaleCallLocaliser.cachedTable(L.ssp, L.D, -Hs(i,3), zl(l), L.o, rmax);
                    TB.T{i,l} = tb.T1; TB.dr = tb.dr;
                end
            end
        end

        function T = locTimes(obj, C, Hs, TB)
            % Localiser travel-time model: K candidates x n hydrophones
            if isempty(TB)
                T = OceanEnvironment.modelTravelTime(obj.cfg.locSSP, C, Hs);
                return;
            end
            % ray time = exact straight-ray time + interpolated ray residual
            K = size(C,1); n = size(Hs,1);
            T = OceanEnvironment.modelTravelTime(obj.cfg.locRay.ssp, C, Hs);
            nl = numel(TB.zl);
            for i = 1:n
                r = sqrt((C(:,1) - Hs(i,1)).^2 + (C(:,2) - Hs(i,2)).^2);
                if nl == 1
                    dT = WhaleCallLocaliser.tableLookup(TB.T{i,1}, TB.dr, r);
                else
                    zs = -C(:,3);
                    u = (zs - TB.zl(1)) / (TB.zl(2) - TB.zl(1));
                    l0 = min(max(floor(u) + 1, 1), nl - 1);
                    w = min(max(u - (l0 - 1), 0), 1);
                    Ta = nan(K,1); Tb = nan(K,1);
                    for l = unique(l0)'
                        m = l0 == l;
                        Ta(m) = WhaleCallLocaliser.tableLookup(TB.T{i,l}, TB.dr, r(m));
                        Tb(m) = WhaleCallLocaliser.tableLookup(TB.T{i,l+1}, TB.dr, r(m));
                    end
                    dT = (1 - w).*Ta + w.*Tb;
                    one = ~isfinite(Ta) & isfinite(Tb); dT(one) = Tb(one);
                    one = isfinite(Ta) & ~isfinite(Tb); dT(one) = Ta(one);
                end
                dT(~isfinite(dT)) = 0;       % shadow zone in the table: straight ray
                T(:,i) = T(:,i) + dT;
            end
        end
    end

    methods (Static)
        function S = cachedPaths(ssp, D, zr, zs, o)
            zr = round(zr*1000)/1000; zs = round(zs*1000)/1000;
            sig = OceanRayModel.sspSignature(ssp);
            ok = sprintf('%d_%.3g_%.3g_%d_%g_%g_%g', o.nRays, o.maxAngle, o.BL, o.maxBottom, o.rmax, o.causticCapDb, o.nLayers);
            fkey = sprintf('F|%s|%.2f|%.3f|%s', sig, D, zr, ok);
            F = wahi_ray_cache(fkey, @() OceanRayModel.fan(ssp, D, zr, o));
            S = wahi_ray_cache(sprintf('P|%s|%.3f', fkey, zs), @() OceanRayModel.paths(F, zs, o));
        end

        function tb = cachedTable(ssp, D, zr, zs, o, rmax)
            zr = round(zr*2)/2; zs = round(zs*10)/10;
            o.rmax = max(o.rmax, rmax);
            o.rmax = 10000*ceil(o.rmax/10000);
            sig = OceanRayModel.sspSignature(ssp);
            key = sprintf('T|%s|%.2f|%.1f|%.1f|%d_%g_%g', sig, D, zr, zs, o.nRays, o.rmax, o.BL);
            tb = wahi_ray_cache(key, @() WhaleCallLocaliser.buildTable(ssp, D, zr, zs, o));
        end

        function tb = buildTable(ssp, D, zr, zs, o)
            % First-arrival table stored as a RESIDUAL relative to the
            % straight-ray time (smooth, so interpolation in range and depth
            % is accurate); NaN = shadow zone
            sig = OceanRayModel.sspSignature(ssp);
            ok = sprintf('%d_%.3g_%.3g_%d_%g_%g_%g', o.nRays, o.maxAngle, o.BL, o.maxBottom, o.rmax, o.causticCapDb, o.nLayers);
            F = wahi_ray_cache(sprintf('F|%s|%.2f|%.3f|%s', sig, D, zr, ok), @() OceanRayModel.fan(ssp, D, zr, o));
            S = OceanRayModel.paths(F, zs, o);
            [T1, dr] = OceanRayModel.firstArrivalTable(S, o.rmax);
            r = (0:numel(T1)-1)' * dr;
            Ts = OceanEnvironment.modelTravelTime(ssp, [r zeros(size(r)) -zs*ones(size(r))], [0 0 -zr]);
            tb = struct('T1', T1 - Ts, 'dr', dr);
        end

        function t = tableLookup(T1, dr, r)
            n = numel(T1);
            u = r(:)/dr;
            i0 = floor(u) + 1;
            fr = u - (i0 - 1);
            out = i0 >= n;
            i0 = min(i0, n - 1);
            t = T1(i0).*(1 - fr) + T1(i0 + 1).*fr;
            t(out) = NaN;
        end

        function cfg = fillDefaults(cfg)
            D = struct('snrMode','sonar','snrFixed',20,'SL',180,'NL',90,'DT',10,'minDetect',4, ...
                'timingSigma',0,'posSigma',0,'solveMode','2d','assumedDepth',50, ...
                'tdoaMethod','waveform','searchRmax',100e3,'seed',NaN,'useTrueDepth',false, ...
                'propModel','straight','locModel','straight');
            fn = fieldnames(D);
            for k = 1:numel(fn)
                if ~isfield(cfg, fn{k}) || isempty(cfg.(fn{k})), cfg.(fn{k}) = D.(fn{k}); end
            end
            if ~isfield(cfg,'locSSP') || isempty(cfg.locSSP)
                cfg.locSSP = OceanEnvironment.constantSSP(1500);
            end
            if ~isfield(cfg,'band') || isempty(cfg.band)
                [a, b] = WhaleCallLibrary.estimateBand(cfg.call, cfg.fs);
                cfg.band = [a b];
            end
            if strcmpi(cfg.propModel,'ray')
                if ~isfield(cfg,'ray') || isempty(cfg.ray), error('propModel ''ray'' needs cfg.ray'); end
                cfg.ray.o = OceanRayModel.fillOpts(cfg.ray.o);
                if ~isfield(cfg.ray,'arrivals'), cfg.ray.arrivals = 'all'; end
                if ~isfield(cfg.ray,'window'), cfg.ray.window = 1; end
            end
            if strcmpi(cfg.locModel,'ray')
                if ~isfield(cfg,'locRay') || isempty(cfg.locRay), error('locModel ''ray'' needs cfg.locRay'); end
                cfg.locRay.o = OceanRayModel.fillOpts(cfg.locRay.o);
                if ~isfield(cfg.locRay,'nZ'), cfg.locRay.nZ = 12; end
            end
            need = 3; if strcmpi(cfg.solveMode,'3d'), need = 4; end
            cfg.minDetect = max(cfg.minDetect, need);
        end

        function Y = envelopeSpectra(Y)
            N = size(Y,1);
            h = zeros(N,1); h(1) = 1; h(2:N/2) = 2; h(N/2+1) = 1;
            for j = 1:size(Y,2)
                e = abs(ifft(Y(:,j).*h));
                e = e - mean(e);
                Y(:,j) = fft(e);
            end
        end

        function out = emptyResult()
            out = struct('flag','', 'est',NaN(1,3), 'err3',NaN, 'errH',NaN, 'errZ',NaN, ...
                'errBrg',NaN, 'errRange',NaN, 'nDet',0, 'minSNR',NaN, 'maxSNR',NaN);
        end

        function out = fillErrors(out, est, src, ctr)
            out.est = est;
            out.err3 = norm(est - src);
            out.errH = norm(est(1:2) - src(1:2));
            out.errZ = abs(est(3) - src(3));
            bT = atan2d(src(1)-ctr(1), src(2)-ctr(2));
            bE = atan2d(est(1)-ctr(1), est(2)-ctr(2));
            out.errBrg = abs(mod(bE - bT + 180, 360) - 180);
            out.errRange = abs(norm(est(1:2)-ctr(1:2)) - norm(src(1:2)-ctr(1:2)));
        end
    end
end
