% OceanRayModel - range-independent acoustic ray tracer for WAH-i Ocean.
%
% Adds refraction, ducting / convergence zones, shadow zones, and surface and
% bottom multipath. The water column 0..D is split into thin layers, each with
% a LINEAR sound-speed gradient, so every ray segment is an exact circular arc.
%
%   Snell's law (ray parameter)   p = cos(theta) / c(z)     (theta = grazing angle)
%   layer of thickness h, gradient g = dc/dz, angles theta1 (top), theta2 (bottom):
%     horizontal range            dr = (sin th1 - sin th2) / (p g)
%     travel time                 dt = (1/g) ln[ (1 + sin th1) cos th2 / ((1 + sin th2) cos th1) ]
%     g -> 0 (iso-speed layer)    dr = h cot th,   dt = h / (c sin th)
%   turning point                 inside the layer where p c(z) = 1
%   surface                       pressure-release reflection, coefficient -1
%   bottom                        specular reflection, loss BL dB per bounce
%
% Because the ocean is horizontally stratified, every ray is PERIODIC: after
% one cycle (down to its lower turning point / bottom, up to its upper turning
% point / surface, and back) it repeats with cycle range Rc and cycle time Tc.
% Each ray is therefore integrated once over depth; crossings of the source
% depth at any range are r = r0 + n Rc, t = t0 + n Tc (n = 0, 1, 2, ...).
%
% Rays are traced from the HYDROPHONE depth (reciprocity). Eigenrays to a
% source at horizontal range R are found by interpolating between neighbouring
% launch angles of the same path family (launch direction, crossing direction,
% cycle number, bounce counts). Intensity from ray-tube spreading:
%
%   TL_spread = 10 log10( R |dr/dtheta0| sin(theta_s) / cos(theta0) )   (dB re 1 m)
%
% which reduces to 20 log10(L) for straight rays. Caustics (dr/dtheta0 -> 0)
% are capped at causticCapDb above spherical spreading. Each eigenray also
% gets n_bottom * BL and seawater absorption alpha(f) * pathLength. Total
% received level = incoherent sum over eigenrays.
%
% Limitations: range-independent (one sound-speed profile and one flat bottom
% for the whole field, normally taken at the array), geometric acoustics only
% (no diffraction into shadow zones, no beam displacement, no bottom
% penetration / sediment model, no frequency dependence except absorption).

classdef OceanRayModel
    methods (Static)

        function o = defaults()
            o = struct('nRays', 400, 'maxAngle', 89.95, 'minAngle', 1e-3, 'maxBottom', 12, ...
                'maxSurface', 400, 'maxCycles', 300, 'BL', 6, 'causticCapDb', 20, ...
                'nLayers', 200, 'rmax', 120e3);
        end

        function o = fillOpts(o)
            d = OceanRayModel.defaults();
            fn = fieldnames(d);
            for k = 1:numel(fn)
                if ~isfield(o, fn{k}) || isempty(o.(fn{k})), o.(fn{k}) = d.(fn{k}); end
            end
        end

        %% ------------------------------------------------------------------
        function F = fan(ssp, D, zr, o)
            % Ray fan launched from hydrophone depth zr (m) in water of depth D.
            o = OceanRayModel.fillOpts(o);
            D = max(D, 1);
            zr = min(max(zr, 0.1), D - 0.1);
            dz = max(D / o.nLayers, 0.5);
            zg = unique([(0:dz:D)'; D; zr]);
            cg = OceanEnvironment.uniformLookup(ssp.Sz, ssp.cz, zg);
            nz = numel(zg);
            ir = find(abs(zg - zr) < 1e-9, 1);
            cr = cg(ir);

            % 30 % of the rays log-spaced near horizontal (long-range direct and
            % ducted paths), the rest uniform up to near-vertical
            nLog = max(10, round(0.3*o.nRays));
            thd = unique([logspace(log10(o.minAngle), log10(3), nLog), linspace(3, o.maxAngle, o.nRays - nLog + 1)]);
            th = thd(:) * pi/180;
            p = cos(th) / cr;
            nR = numel(p);
            PC = p * cg';                         % nR x nz
            closed = PC >= 1 - 1e-12;

            % vertical extent around zr
            kk = repmat(1:nz, nR, 1);
            cAbove = closed(:, 1:ir) .* kk(:, 1:ir);
            lastClosedAbove = max(cAbove, [], 2);
            kTop = lastClosedAbove + 1;
            turnTop = lastClosedAbove > 0;
            cBelow = kk(:, ir:nz); cBelow(~closed(:, ir:nz)) = inf;
            firstClosedBelow = min(cBelow, [], 2);
            turnBot = isfinite(firstClosedBelow);
            kBot = nz*ones(nR,1); kBot(turnBot) = firstClosedBelow(turnBot) - 1;

            % full layers (open at both ends)
            h = diff(zg)';                           % 1 x (nz-1)
            ca = cg(1:end-1)'; cb = cg(2:end)';
            [LR, LT] = OceanRayModel.seg(repmat(p,1,nz-1), repmat(ca,nR,1), repmat(cb,nR,1), repmat(h,nR,1));
            both = ~closed(:,1:end-1) & ~closed(:,2:end);
            LR(~both) = 0; LT(~both) = 0;
            CR = [zeros(nR,1) cumsum(LR, 2)];
            CT = [zeros(nR,1) cumsum(LT, 2)];

            % partial turning layers
            PT = zeros(nR,1); PTt = zeros(nR,1); zTurnTop = zeros(nR,1);
            PB = zeros(nR,1); PBt = zeros(nR,1); zTurnBot = D*ones(nR,1);
            j = find(turnTop);
            if ~isempty(j)
                k1 = kTop(j) - 1;
                g = (cg(k1+1) - cg(k1)) ./ (zg(k1+1) - zg(k1));
                zt = zg(k1) + (1./p(j) - cg(k1)) ./ g;
                zTurnTop(j) = zt;
                [PT(j), PTt(j)] = OceanRayModel.seg(p(j), 1./p(j), cg(k1+1), zg(k1+1) - zt);
            end
            j = find(turnBot);
            if ~isempty(j)
                k1 = kBot(j);
                g = (cg(k1+1) - cg(k1)) ./ (zg(k1+1) - zg(k1));
                zt = zg(k1) + (1./p(j) - cg(k1)) ./ g;
                zTurnBot(j) = zt;
                [PB(j), PBt(j)] = OceanRayModel.seg(p(j), cg(k1), 1./p(j), zt - zg(k1));
            end

            rowi = (1:nR)';
            CRtop = CR(sub2ind(size(CR), rowi, kTop)); CTtop = CT(sub2ind(size(CT), rowi, kTop));
            Rbot = PT + CR(sub2ind(size(CR), rowi, kBot)) - CRtop + PB;
            Tbot = PTt + CT(sub2ind(size(CT), rowi, kBot)) - CTtop + PBt;

            F = struct();
            F.zr = zr; F.D = D; F.zg = zg; F.cg = cg; F.ir = ir; F.cr = cr;
            F.th0 = th; F.p = p; F.kTop = kTop; F.kBot = kBot;
            F.turnTop = turnTop; F.turnBot = turnBot; F.zTurnTop = zTurnTop; F.zTurnBot = zTurnBot;
            F.CR = CR; F.CT = CT; F.CRtop = CRtop; F.CTtop = CTtop;
            F.PT = PT; F.PTt = PTt;
            F.Rc = 2*Rbot; F.Tc = 2*Tbot;
            F.RdR = PT + CR(sub2ind(size(CR), rowi, ir*ones(nR,1))) - CRtop;
            F.TdR = PTt + CT(sub2ind(size(CT), rowi, ir*ones(nR,1))) - CTtop;
            F.o = o;
        end

        function [dr, dt] = seg(p, ca, cb, h)
            % Exact range/time for a linear-gradient layer from c = ca to c = cb
            % (thickness h > 0). An end with c = 1/p is a turning point.
            dr = zeros(size(h)); dt = zeros(size(h));
            ok = h > 1e-9;
            if ~any(ok(:)), return; end
            p = p(ok); ca = ca(ok); cb = cb(ok); hh = h(ok);
            g = (cb - ca) ./ hh;
            xa = min(p.*ca, 1); xb = min(p.*cb, 1);
            sa = sqrt(max(1 - xa.^2, 0)); sb = sqrt(max(1 - xb.^2, 0));
            r = (sa - sb) ./ (p.*g);
            t = (log((1 + sa)./xa) - log((1 + sb)./xb)) ./ g;
            sm = abs(g) < 1e-7;
            if any(sm)
                s = max(0.5*(sa(sm) + sb(sm)), 1e-12);
                x = 0.5*(xa(sm) + xb(sm));
                r(sm) = hh(sm) .* x ./ s;
                t(sm) = hh(sm) ./ (0.5*(ca(sm) + cb(sm)) .* s);
            end
            dr(ok) = r; dt(ok) = t;
        end

        function [Rd, Td, ss, valid] = RdAt(F, zs)
            % Canonical down-leg range/time from the upper turning point (or
            % surface) to depth zs, the grazing-angle sine there, and validity.
            zg = F.zg; cg = F.cg; p = F.p; nR = numel(p);
            zs = min(max(zs, 0), F.D);
            k = find(zg <= zs + 1e-12, 1, 'last');
            if k >= numel(zg), k = numel(zg) - 1; end
            cs = cg(k) + (cg(k+1) - cg(k)) * (zs - zg(k)) / (zg(k+1) - zg(k));
            ss = sqrt(max(1 - (p*cs).^2, 0));
            zTop = F.zTurnTop; zTop(~F.turnTop) = 0;
            zBot = F.zTurnBot; zBot(~F.turnBot) = F.D;
            valid = (p*cs < 1 - 1e-12) & zs >= zTop - 1e-9 & zs <= zBot + 1e-9;
            Rd = nan(nR,1); Td = nan(nR,1);
            % partial top layer (zs above the first open level)
            a = valid & zs < zg(F.kTop) - 1e-12;
            if any(a)
                [Rd(a), Td(a)] = OceanRayModel.seg(p(a), 1./p(a), cs*ones(sum(a),1), zs - F.zTurnTop(a));
            end
            b = valid & ~a;
            if any(b)
                kb = min(k, F.kBot(b));
                rowi = find(b);
                base = F.PT(b) + F.CR(sub2ind(size(F.CR), rowi, kb)) - F.CRtop(b);
                baseT = F.PTt(b) + F.CT(sub2ind(size(F.CT), rowi, kb)) - F.CTtop(b);
                [dr, dt] = OceanRayModel.seg(p(b), cg(kb), cs*ones(sum(b),1), zs - zg(kb));
                Rd(b) = base + dr; Td(b) = baseT + dt;
            end
            valid = valid & isfinite(Rd) & F.Rc > 1e-6;
        end

        %% ------------------------------------------------------------------
        function S = paths(F, zs, o)
            % All eigenray segments from hydrophone depth F.zr to source depth zs
            o = OceanRayModel.fillOpts(o);
            [RdS, TdS, sS, valid] = OceanRayModel.RdAt(F, zs);
            Rc = F.Rc; Tc = F.Tc;
            surfFlag = ~F.turnTop; botFlag = ~F.turnBot;
            rLim = 1.15*o.rmax + 1000;
            minRc = min(Rc(valid & Rc > 0));
            if isempty(minRc), minRc = rLim; end
            nCyc = min(o.maxCycles, ceil(rLim/minRc) + 1);
            nvec = 0:nCyc-1;
            S = struct('rA',[],'rB',[],'tA',[],'tB',[],'thA',[],'thB',[],'sA',[],'sB',[],'nb',[],'ns',[]);
            for d = [1 -1]
                if d > 0, s0 = F.RdR; t0 = F.TdR; else, s0 = Rc - F.RdR; t0 = Tc - F.TdR; end
                for x = [1 -1]
                    if x > 0, cr = RdS; ct = TdS; else, cr = Rc - RdS; ct = Tc - TdS; end
                    dfr = cr - s0; dft = ct - t0;
                    w = dfr < -1e-9;
                    offr = dfr + w.*Rc; offt = dft + w.*Tc;
                    r = offr + Rc*nvec; t = offt + Tc*nvec;          % nR x nCyc
                    e = s0 + r;
                    nBot = floor((e - Rc/2)./Rc + 1e-9) - floor((s0 - Rc/2)./Rc + 1e-9);
                    nTop = floor(e./Rc + 1e-9) - floor(s0./Rc + 1e-9);
                    nb = nBot .* botFlag; ns = nTop .* surfFlag;
                    ok = valid & r <= rLim & nb <= o.maxBottom & ns <= o.maxSurface & r >= 0;
                    % adjacent rays of the same family
                    pair = ok(1:end-1,:) & ok(2:end,:) & nb(1:end-1,:) == nb(2:end,:) & ns(1:end-1,:) == ns(2:end,:) ...
                        & (surfFlag(1:end-1) == surfFlag(2:end)) & (botFlag(1:end-1) == botFlag(2:end));
                    [jj, nn] = find(pair);
                    if isempty(jj), continue; end
                    iA = sub2ind(size(r), jj, nn); iB = sub2ind(size(r), jj+1, nn);
                    S.rA = [S.rA; r(iA)]; S.rB = [S.rB; r(iB)];
                    S.tA = [S.tA; t(iA)]; S.tB = [S.tB; t(iB)];
                    S.thA = [S.thA; F.th0(jj)]; S.thB = [S.thB; F.th0(jj+1)];
                    S.sA = [S.sA; sS(jj)]; S.sB = [S.sB; sS(jj+1)];
                    S.nb = [S.nb; nb(iA)]; S.ns = [S.ns; ns(iA)];
                end
            end
            keep = S.rA ~= S.rB;
            fn = fieldnames(S);
            for k = 1:numel(fn), S.(fn{k}) = S.(fn{k})(keep); end
            S.zr = F.zr; S.zs = zs; S.o = o;
            % straight-path mean slowness between zr and zs (direct-path fallback)
            zq = linspace(min(F.zr, zs), max(F.zr, zs), 200)';
            cq = interp1(F.zg, F.cg, zq, 'linear', 'extrap');
            S.slowDirect = mean(1 ./ cq);
            hasDirect = S.nb == 0 & S.ns == 0;
            S.directMaxR = 0;
            if any(hasDirect), S.directMaxR = max(max(S.rA(hasDirect)), max(S.rB(hasDirect))); end
            S.directMinR = inf;
            if any(hasDirect), S.directMinR = min(min(S.rA(hasDirect)), min(S.rB(hasDirect))); end
        end

        function A = eigenrays(S, R, f_Hz, absorption)
            % Eigenrays at horizontal range R. Returns struct with t (s),
            % TL (dB re 1 m), sgn (+/-1), nb, ns, th0 (launch angle at the hydrophone)
            if nargin < 4, absorption = true; end
            o = S.o;
            m = (R - S.rA).*(R - S.rB) <= 0;
            A = struct('t',[],'TL',[],'sgn',[],'nb',[],'ns',[],'th0',[]);
            L = sqrt(R^2 + (S.zs - S.zr)^2);
            % Direct-path fallback where the fan cannot resolve it: nearly
            % vertical (R below the steepest direct ray) or source and receiver
            % at (almost) the same depth (grazing ray degenerates)
            needDirect = ~any(m & S.nb == 0 & S.ns == 0) && ...
                (R < S.directMinR || abs(S.zs - S.zr) < max(0.5, 1e-3*R));
            if ~any(m) && ~needDirect, return; end
            if needDirect
                tD = L*S.slowDirect;
                TLD = 20*log10(max(L,1));
                if absorption, TLD = TLD + OceanEnvironment.thorp(f_Hz/1000)*L/1000; end
            end
            if ~any(m)
                A.t = tD; A.TL = TLD; A.sgn = 1; A.nb = 0; A.ns = 0; A.th0 = atan2(abs(S.zs - S.zr), max(R,1e-6));
                return;
            end
            rA = S.rA(m); rB = S.rB(m);
            u = (R - rA)./(rB - rA);
            t = S.tA(m) + u.*(S.tB(m) - S.tA(m));
            th = S.thA(m) + u.*(S.thB(m) - S.thA(m));
            sS = max(S.sA(m) + u.*(S.sB(m) - S.sA(m)), 1e-4);
            drdth = abs(rB - rA)./max(abs(S.thB(m) - S.thA(m)), 1e-12);
            Rm = max(R, 1);
            TL = 10*log10(Rm .* drdth .* sS ./ max(cos(th), 1e-6));
            TL = max(TL, 20*log10(max(L,1)) - o.causticCapDb);
            TL = TL + S.nb(m)*o.BL;
            if absorption
                TL = TL + OceanEnvironment.thorp(f_Hz/1000) .* (t*1500) / 1000;
            end
            A.t = t; A.TL = TL; A.sgn = (-1).^S.ns(m); A.nb = S.nb(m); A.ns = S.ns(m); A.th0 = th;
            if needDirect
                A.t = [tD; A.t]; A.TL = [TLD; A.TL]; A.sgn = [1; A.sgn]; A.nb = [0; A.nb]; A.ns = [0; A.ns];
                A.th0 = [atan2(abs(S.zs - S.zr), max(R,1e-6)); A.th0];
            end
            [A.t, ord] = sort(A.t);
            A.TL = A.TL(ord); A.sgn = A.sgn(ord); A.nb = A.nb(ord); A.ns = A.ns(ord); A.th0 = A.th0(ord);
        end

        function [T1, dr] = firstArrivalTable(S, rmax, dr)
            % First-arrival time on a uniform range grid r = (0:n-1)*dr (NaN = shadow)
            if nargin < 3 || isempty(dr), dr = max(rmax/4000, 5); end
            n = floor(rmax/dr) + 1;
            T1 = nan(n,1);
            if isempty(S.rA), return; end
            lo = min(S.rA, S.rB); hi = max(S.rA, S.rB);
            i1 = max(ceil(lo/dr) + 1, 1); i2 = min(floor(hi/dr) + 1, n);
            cnt = max(i2 - i1 + 1, 0);
            tot = sum(cnt);
            if tot == 0, return; end
            sIdx = repelem((1:numel(cnt))', cnt);
            first = cumsum(cnt) - cnt;
            gIdx = i1(sIdx) + ((1:tot)' - 1 - first(sIdx));
            rr = (gIdx - 1)*dr;
            u = (rr - S.rA(sIdx)) ./ (S.rB(sIdx) - S.rA(sIdx));
            tt = S.tA(sIdx) + u .* (S.tB(sIdx) - S.tA(sIdx));
            T1 = accumarray(gIdx, tt, [n 1], @min, NaN);
        end

        %% ------------------------------------------------------------------
        function P = rayPaths(F, anglesDeg, rmax)
            % Polylines [r z] for plotting rays launched up and down at the
            % requested grazing angles from the hydrophone depth.
            P = {};
            for a = anglesDeg(:)'
                [~, j] = min(abs(F.th0*180/pi - abs(a)));
                Rc = F.Rc(j);
                if ~(Rc > 1e-6), continue; end
                ks = F.kTop(j):F.kBot(j);
                rd = F.PT(j) + F.CR(j, ks) - F.CRtop(j);
                zTop = 0; if F.turnTop(j), zTop = F.zTurnTop(j); end
                zBot = F.D; if F.turnBot(j), zBot = F.zTurnBot(j); end
                cyc = [0 rd Rc/2 Rc - fliplr(rd) Rc; zTop F.zg(ks)' zBot fliplr(F.zg(ks)') zTop];
                for d = [1 -1]
                    if d > 0, s0 = F.RdR(j); else, s0 = Rc - F.RdR(j); end
                    nC = ceil((rmax + s0)/Rc) + 1;
                    nC = min(nC, 2000);
                    rr = []; zz = [];
                    for n = 0:nC-1
                        rr = [rr cyc(1,:) + n*Rc]; %#ok<AGROW>
                        zz = [zz cyc(2,:)]; %#ok<AGROW>
                    end
                    rr = rr - s0;
                    keep = rr > 0 & rr <= rmax;
                    P{end+1} = [[0 rr(keep)]' [F.zr zz(keep)]']; %#ok<AGROW>
                end
            end
        end

        function sig = sspSignature(ssp)
            c = ssp.cz(:);
            sig = sprintf('%s_%d_%.8g_%.8g', ssp.type, numel(c), sum(c), sum(c .* (1:numel(c))') / numel(c));
        end
    end
end
