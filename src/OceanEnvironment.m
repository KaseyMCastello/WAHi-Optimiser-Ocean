% OceanEnvironment - bathymetry, sound speed and propagation helpers for WAH-i Ocean.
%
% Coordinate convention used throughout WAH-i Ocean
% -------------------------------------------------
%   x : metres East of the local origin
%   y : metres North of the local origin
%   z : metres ELEVATION relative to the sea surface (z = 0 at the surface,
%       negative below it). Depth = -z.
% Geographic inputs (lon/lat) are converted to local x/y with an
% equirectangular projection about (lon0, lat0), which is accurate to well
% under 0.1 % over the 100 km scales used here.
%
% An environment is a plain struct (see OceanEnvironment.default) so it can be
% saved in sessions and passed to parfor workers:
%   env.flatDepth    water depth (m) used when no bathymetry is loaded
%   env.bathy        [] or struct(x, y, depth, isGeo, lon0, lat0, file, F)
%   env.lat0/lon0    geographic origin (NaN if unknown)
%   env.ssp          struct(type = 'constant'|'profile'|'field', c0, ...)
%   env.prop         struct(spreading, transR, absorption)
%
% All methods are static. Main entry points:
%   env  = OceanEnvironment.default()
%   B    = OceanEnvironment.loadBathymetry(file, lon0, lat0)
%   d    = OceanEnvironment.waterDepth(env, x, y)
%   prof = OceanEnvironment.loadProfile(file)
%   fld  = OceanEnvironment.loadNCOM(file, lon0, lat0, tIdx)
%   t    = OceanEnvironment.travelTime(env, src, rcv)
%   TL   = OceanEnvironment.transmissionLoss(env, r, f)
%
% Sound propagation model: straight-line rays. For a depth-dependent profile
% the travel time along a straight path is exact for that path
% (t = L * mean slowness over the depth interval). Refraction, ducting,
% shadow zones and multipath are NOT modelled; see README for caveats.

classdef OceanEnvironment
    methods (Static)

        %% ------------------------------------------------------------------
        function env = default()
            env = struct();
            env.flatDepth = 1000;
            env.bathy = [];
            env.lat0 = NaN;
            env.lon0 = NaN;
            env.ssp = OceanEnvironment.constantSSP(1500);
            env.prop = struct('spreading','hybrid','transR',1000,'absorption',true);
        end

        function ssp = constantSSP(c0)
            ssp = struct('type','constant','c0',c0,'z',[0;12000],'c',[c0;c0], ...
                'Sz',[],'S',[],'field',[],'file','');
            ssp = OceanEnvironment.prepareProfile(ssp);
        end

        %% ------------------------------------------------------------------
        %  Geographic helpers
        function [x, y] = geo2local(lon, lat, lon0, lat0)
            x = (lon - lon0) .* cosd(lat0) * 111320;
            y = (lat - lat0) * 110574;
        end

        function [lon, lat] = local2geo(x, y, lon0, lat0)
            lon = lon0 + x ./ (cosd(lat0) * 111320);
            lat = lat0 + y / 110574;
        end

        %% ------------------------------------------------------------------
        %  Bathymetry
        function B = loadBathymetry(file, lon0, lat0)
            % Loads a bathymetry grid or scattered XYZ file.
            % Supported: .csv/.txt/.xyz/.dat (3 columns x y z or lon lat z),
            %            .mat (X,Y,Z or lon,lat,z/depth/elevation or Nx3 matrix),
            %            .nc/.grd/.nc4 (GEBCO/ETOPO/SRTM style lon/lat/elevation).
            % The z value may be depth (positive down) or elevation (negative
            % below sea level): if the median is negative it is treated as
            % elevation and flipped. Geographic coordinates are detected
            % automatically and converted to local metres about (lon0, lat0);
            % if lon0/lat0 are NaN the grid centre is used.
            if nargin < 2, lon0 = NaN; end
            if nargin < 3, lat0 = NaN; end
            if ~exist(file,'file'), error('Bathymetry file not found: %s', file); end
            [~,~,ext] = fileparts(file);
            ext = lower(ext);
            isGeoHint = [];
            switch ext
                case {'.nc','.grd','.nc4','.cdf'}
                    [gx, gy, gz, isGeoHint] = OceanEnvironment.readNetcdfGrid(file);
                case '.mat'
                    [gx, gy, gz, isGeoHint] = OceanEnvironment.readMatGrid(file);
                otherwise
                    A = OceanEnvironment.readNumericTable(file);
                    if size(A,2) < 3
                        error('Bathymetry text file needs 3 columns (x y z or lon lat z).');
                    end
                    [gx, gy, gz] = OceanEnvironment.xyzToGrid(A(:,1), A(:,2), A(:,3));
            end

            gx = double(gx(:))'; gy = double(gy(:))'; gz = double(gz);
            % ensure increasing axes
            if numel(gx) > 1 && gx(2) < gx(1), gx = fliplr(gx); gz = fliplr(gz); end
            if numel(gy) > 1 && gy(2) < gy(1), gy = fliplr(gy); gz = flipud(gz); end

            % geographic?
            if isempty(isGeoHint)
                isGeo = all(abs(gy) <= 90) && all(gx >= -180) && all(gx <= 360) && ...
                    (max(gx)-min(gx)) < 30 && (max(gy)-min(gy)) < 30;
            else
                isGeo = isGeoHint;
            end

            % depth sign
            zfin = gz(isfinite(gz));
            flipped = false;
            if ~isempty(zfin) && median(zfin) < 0
                gz = -gz; flipped = true;
            end

            B = struct();
            B.file = file;
            B.isGeo = isGeo;
            B.signFlipped = flipped;
            if isGeo
                gx(gx > 180) = gx(gx > 180) - 360;           % 0..360 -> -180..180
                [gx, ordx] = sort(gx); gz = gz(:, ordx);
                if ~isfinite(lon0), lon0 = mean([min(gx) max(gx)]); end
                if ~isfinite(lat0), lat0 = mean([min(gy) max(gy)]); end
                if lon0 > 180, lon0 = lon0 - 360; end
                B.lonVec = gx; B.latVec = gy;
                [xx, ~] = OceanEnvironment.geo2local(gx, lat0*ones(size(gx)), lon0, lat0);
                [~, yy] = OceanEnvironment.geo2local(lon0*ones(size(gy)), gy, lon0, lat0);
                B.x = xx; B.y = yy;
            else
                B.x = gx; B.y = gy;
                B.lonVec = []; B.latVec = [];
            end
            B.lon0 = lon0; B.lat0 = lat0;
            B.depth = gz;   % ny x nx, positive down (negative = land)
            B = OceanEnvironment.finaliseBathy(B);
        end

        function B = finaliseBathy(B)
            % fill NaNs with nearest valid value so interpolation is safe
            D = B.depth;
            if any(~isfinite(D(:)))
                good = isfinite(D);
                if ~any(good(:)), error('Bathymetry grid contains no finite values.'); end
                [YY, XX] = ndgrid(1:size(D,1), 1:size(D,2));
                try
                    Fn = scatteredInterpolant(XX(good), YY(good), D(good), 'nearest', 'nearest');
                    D(~good) = Fn(XX(~good), YY(~good));
                catch
                    D(~good) = griddata(XX(good), YY(good), D(good), XX(~good), YY(~good), 'nearest');
                end
                B.depth = D;
            end
            B.F = griddedInterpolant({B.y(:), B.x(:)}, B.depth, 'linear', 'nearest');
            B.minDepth = min(B.depth(:));
            B.maxDepth = max(B.depth(:));
        end

        function B = flipBathySign(B)
            B.depth = -B.depth;
            B.signFlipped = ~B.signFlipped;
            B = OceanEnvironment.finaliseBathy(B);
        end

        function d = waterDepth(env, x, y)
            % Water depth (m, positive) at points. Land returns <= 0.
            if isempty(env.bathy)
                d = env.flatDepth * ones(size(x));
            else
                d = env.bathy.F(y, x);
            end
        end

        function dmax = maxWaterDepth(env)
            if isempty(env.bathy), dmax = env.flatDepth;
            else, dmax = max(env.bathy.maxDepth, 1);
            end
        end

        %% ------------------------------------------------------------------
        %  Sound speed
        function c = mackenzie(T, S, D)
            % Mackenzie (1981) nine-term equation. T degC, S psu, D depth m.
            c = 1448.96 + 4.591*T - 5.304e-2*T.^2 + 2.374e-4*T.^3 + ...
                1.340*(S-35) + 1.630e-2*D + 1.675e-7*D.^2 - ...
                1.025e-2*T.*(S-35) - 7.139e-13*T.*D.^3;
        end

        function ssp = loadProfile(file)
            % 1-D sound speed profile from CSV/TXT/MAT.
            %   2 columns: depth (m) and sound speed (m/s), either order
            %   3 columns: depth (m), temperature (degC), salinity (psu)
            [~,~,ext] = fileparts(file);
            if strcmpi(ext,'.mat')
                Sm = load(file);
                fn = fieldnames(Sm);
                A = [];
                for k = 1:numel(fn)
                    v = Sm.(fn{k});
                    if isnumeric(v) && size(v,2) >= 2 && size(v,1) >= 2, A = double(v); break; end
                end
                if isempty(A), error('No N x 2 numeric matrix found in %s', file); end
            else
                A = OceanEnvironment.readNumericTable(file);
            end
            A = A(all(isfinite(A(:,1:min(3,size(A,2)))),2),:);
            if size(A,1) < 2, error('Sound speed profile needs at least 2 rows.'); end
            if size(A,2) >= 3 && all(A(:,2) > -5 & A(:,2) < 40) && all(A(:,3) >= 0 & A(:,3) < 45) ...
                    && ~all(A(:,2) > 1300)
                z = abs(A(:,1));
                c = OceanEnvironment.mackenzie(A(:,2), A(:,3), z);
                note = 'depth/T/S (Mackenzie 1981)';
            else
                isC1 = mean(A(:,1) > 1300 & A(:,1) < 1700);
                isC2 = mean(A(:,2) > 1300 & A(:,2) < 1700);
                if isC2 >= isC1
                    z = abs(A(:,1)); c = A(:,2);
                else
                    z = abs(A(:,2)); c = A(:,1);
                end
                note = 'depth/c';
            end
            [z, ia] = unique(z);
            c = c(ia);
            ssp = struct('type','profile','c0',mean(c),'z',z(:),'c',c(:), ...
                'Sz',[],'S',[],'field',[],'file',file,'note',note);
            ssp = OceanEnvironment.prepareProfile(ssp);
        end

        function ssp = prepareProfile(ssp)
            % Precompute cumulative slowness S(d) = int_0^d dz / c(z) on a 1 m grid
            z = ssp.z(:); c = ssp.c(:);
            if z(1) > 0, z = [0; z]; c = [c(1); c]; end
            if z(end) < 12000, z = [z; 12000]; c = [c; c(end)]; end
            dz = 1;
            if z(end) > 6000, dz = 2; end
            zz = (0:dz:z(end))';
            cc = interp1(z, c, zz, 'linear', 'extrap');
            s = 1 ./ cc;
            S = [0; cumsum(0.5*(s(1:end-1)+s(2:end))*dz)];
            ssp.Sz = zz;
            ssp.S = S;
            ssp.cz = cc;
            ssp.c0 = mean(c(z <= max(ssp.z)));
        end

        function env = rebuild(env)
            % Re-create interpolant objects (e.g. after loading from a MAT file)
            if ~isempty(env.bathy), env.bathy = OceanEnvironment.finaliseBathy(env.bathy); end
            if isfield(env,'ssp') && strcmp(env.ssp.type,'field') && ~isempty(env.ssp.field)
                env.ssp.field = OceanEnvironment.rebuildField(env.ssp.field);
            end
        end

        function fld = rebuildField(fld)
            fld.F = griddedInterpolant({fld.x, fld.y, fld.z}, fld.C, 'linear', 'nearest');
        end

        function fld = loadNCOM(file, lon0, lat0, tIdx)
            % Loads a 3-D T/S (or sound speed) field from an NCOM / HYCOM /
            % CMEMS style NetCDF file and converts it to sound speed on a local
            % x/y/depth grid. Variable names are detected automatically.
            if nargin < 2, lon0 = NaN; end
            if nargin < 3, lat0 = NaN; end
            if nargin < 4 || isempty(tIdx), tIdx = 1; end
            info = ncinfo(file);
            names = {info.Variables.Name};
            vLon = OceanEnvironment.pickName(names, {'lon','longitude','x','nav_lon'});
            vLat = OceanEnvironment.pickName(names, {'lat','latitude','y','nav_lat'});
            vDep = OceanEnvironment.pickName(names, {'depth','z','lev','deptht','level'});
            vC   = OceanEnvironment.pickName(names, {'sound_speed','soundspeed','ssp','c','sound_velocity'});
            vT   = OceanEnvironment.pickName(names, {'water_temp','temperature','temp','thetao','potemp','t','to'});
            vS   = OceanEnvironment.pickName(names, {'salinity','salt','so','sal','s'});
            if isempty(vLon) || isempty(vLat) || isempty(vDep)
                error('Could not find lon/lat/depth variables in %s (found: %s).', file, strjoin(names, ', '));
            end
            lon = double(ncread(file, vLon)); lat = double(ncread(file, vLat));
            dep = abs(double(ncread(file, vDep)));
            if ~isvector(lon), lon = lon(:,1); end
            if ~isvector(lat), lat = lat(1,:); end
            lon = lon(:); lat = lat(:); dep = dep(:);

            if ~isempty(vC)
                C = OceanEnvironment.readField4(file, info, vC, {vLon,vLat,vDep}, tIdx);
                note = sprintf('sound speed variable "%s"', vC);
            else
                if isempty(vT) || isempty(vS)
                    error('Need water temperature and salinity (or sound speed) in %s.', file);
                end
                T = OceanEnvironment.readField4(file, info, vT, {vLon,vLat,vDep}, tIdx);
                Sa = OceanEnvironment.readField4(file, info, vS, {vLon,vLat,vDep}, tIdx);
                if OceanEnvironment.nanmedianLocal(T(:)) > 200, T = T - 273.15; end   % Kelvin
                D3 = repmat(reshape(dep,1,1,[]), numel(lon), numel(lat), 1);
                C = OceanEnvironment.mackenzie(T, Sa, D3);
                note = sprintf('T "%s" / S "%s" (Mackenzie 1981)', vT, vS);
            end
            C(C < 1300 | C > 1700) = NaN;

            % Fill below-seafloor NaNs downward, then fully-empty columns with the mean profile
            meanProf = zeros(numel(dep),1);
            for k = 1:numel(dep)
                v = C(:,:,k); v = v(isfinite(v));
                if isempty(v), meanProf(k) = NaN; else, meanProf(k) = mean(v); end
            end
            for k = 2:numel(dep)
                if ~isfinite(meanProf(k)), meanProf(k) = meanProf(k-1); end
            end
            if ~isfinite(meanProf(1)), error('NCOM field contains no valid sound speeds.'); end
            for i = 1:numel(lon)
                for j = 1:numel(lat)
                    col = squeeze(C(i,j,:));
                    if ~any(isfinite(col))
                        C(i,j,:) = reshape(meanProf,1,1,[]);
                        continue;
                    end
                    for k = 2:numel(col)
                        if ~isfinite(col(k)), col(k) = col(k-1); end
                    end
                    k1 = find(isfinite(col),1);
                    col(1:k1-1) = col(k1);
                    C(i,j,:) = reshape(col,1,1,[]);
                end
            end

            lon(lon > 180) = lon(lon > 180) - 360;          % 0..360 -> -180..180
            if ~isfinite(lon0), lon0 = mean([min(lon) max(lon)]); end
            if ~isfinite(lat0), lat0 = mean([min(lat) max(lat)]); end
            if lon0 > 180, lon0 = lon0 - 360; end
            [lon, il] = sort(lon); [lat, jl] = sort(lat); [dep, kl] = sort(dep);
            C = C(il, jl, kl);
            [x, ~] = OceanEnvironment.geo2local(lon, lat0*ones(size(lon)), lon0, lat0);
            [~, y] = OceanEnvironment.geo2local(lon0*ones(size(lat)), lat, lon0, lat0);

            fld = struct();
            fld.file = file; fld.note = note;
            fld.lon = lon; fld.lat = lat; fld.lon0 = lon0; fld.lat0 = lat0;
            fld.x = x(:); fld.y = y(:); fld.z = dep(:); fld.C = C;
            if numel(x) == 1, fld.x = [x-1e6; x+1e6]; fld.C = cat(1, C, C); end
            if numel(y) == 1, fld.y = [y-1e6; y+1e6]; fld.C = cat(2, fld.C, fld.C); end
            if numel(dep) == 1, fld.z = [dep; dep+1e4]; fld.C = cat(3, fld.C, fld.C); end
            fld.F = griddedInterpolant({fld.x, fld.y, fld.z}, fld.C, 'linear', 'nearest');
            fld.meanProfile = meanProf(kl);
        end

        function ssp = fieldSSP(fld)
            ssp = struct('type','field','c0',mean(fld.meanProfile),'z',fld.z,'c',fld.meanProfile, ...
                'Sz',[],'S',[],'field',fld,'file',fld.file,'note',fld.note);
            ssp = OceanEnvironment.prepareProfile(ssp);
        end

        function prof = profileAt(env, x, y)
            % 1-D profile (struct usable by travel-time model) at a horizontal position
            switch env.ssp.type
                case 'field'
                    fld = env.ssp.field;
                    zq = fld.z(:);
                    cq = fld.F(x*ones(size(zq)), y*ones(size(zq)), zq);
                    prof = struct('type','profile','c0',mean(cq),'z',zq,'c',cq, ...
                        'Sz',[],'S',[],'field',[],'file','','note','profile at array');
                    prof = OceanEnvironment.prepareProfile(prof);
                otherwise
                    prof = env.ssp;
            end
        end

        %% ------------------------------------------------------------------
        %  Travel time
        function t = travelTime(env, src, rcv)
            % Straight-ray travel time (s) from source(s) to receivers.
            %   src : K x 3 or 1 x 3 ; rcv : N x 3  -> t : K x N
            t = OceanEnvironment.modelTravelTime(env.ssp, src, rcv);
        end

        function t = modelTravelTime(ssp, src, rcv)
            % Travel time for a sound-speed model struct (constant / profile / field).
            K = size(src,1); N = size(rcv,1);
            dx = reshape(src(:,1),K,1) - reshape(rcv(:,1),1,N);
            dy = reshape(src(:,2),K,1) - reshape(rcv(:,2),1,N);
            dzv = reshape(src(:,3),K,1) - reshape(rcv(:,3),1,N);
            L = sqrt(dx.^2 + dy.^2 + dzv.^2);
            switch ssp.type
                case 'constant'
                    t = L / ssp.c0;
                case 'profile'
                    ds = max(-src(:,3), 0); dr = max(-rcv(:,3), 0);
                    ds = repmat(reshape(ds,K,1),1,N); dr = repmat(reshape(dr,1,N),K,1);
                    Ss = OceanEnvironment.uniformLookup(ssp.Sz, ssp.S, ds);
                    Sr = OceanEnvironment.uniformLookup(ssp.Sz, ssp.S, dr);
                    dd = dr - ds;
                    slow = (Sr - Ss) ./ dd;
                    small = abs(dd) < 0.5;
                    if any(small(:))
                        cm = OceanEnvironment.uniformLookup(ssp.Sz, ssp.cz, 0.5*(ds(small)+dr(small)));
                        slow(small) = 1 ./ cm;
                    end
                    t = L .* slow;
                case 'field'
                    fld = ssp.field;
                    nS = 24;
                    u = ((1:nS) - 0.5) / nS;
                    t = zeros(K, N);
                    for k = 1:K
                        P0 = src(k,:);
                        D = rcv - P0;  % N x 3
                        XX = P0(1) + D(:,1) * u;   % N x nS
                        YY = P0(2) + D(:,2) * u;
                        ZZ = max(-(P0(3) + D(:,3) * u), 0);
                        cc = fld.F(XX, YY, ZZ);
                        t(k,:) = (L(k,:)' .* mean(1 ./ cc, 2))';
                    end
                otherwise
                    error('Unknown sound speed model type: %s', ssp.type);
            end
        end

        %% ------------------------------------------------------------------
        %  Propagation loss
        function a = thorp(f_kHz)
            % Thorp absorption (dB/km), f in kHz
            f2 = f_kHz.^2;
            a = 0.11*f2./(1+f2) + 44*f2./(4100+f2) + 2.75e-4*f2 + 0.003;
        end

        function TL = transmissionLoss(env, r, f_Hz)
            % One-way transmission loss (dB re 1 m) at slant range r (m).
            r = max(r, 1);
            switch lower(env.prop.spreading)
                case 'spherical',   TL = 20*log10(r);
                case 'cylindrical', TL = 10*log10(r);
                case 'practical',   TL = 15*log10(r);
                otherwise % hybrid: spherical to transR, cylindrical beyond
                    rt = max(env.prop.transR, 1);
                    TL = 20*log10(min(r, rt)) + 10*log10(max(r, rt)/rt);
            end
            if env.prop.absorption
                TL = TL + OceanEnvironment.thorp(f_Hz/1000) .* r / 1000;
            end
        end

        %% ------------------------------------------------------------------
        %  File helpers
        function A = readNumericTable(file)
            % Robust numeric reader for CSV/TXT with optional header/comment lines.
            fid = fopen(file, 'r');
            if fid < 0, error('Cannot open %s', file); end
            cleaner = onCleanup(@() fclose(fid));
            rows = {};
            nc = 0;
            while true
                ln = fgetl(fid);
                if ~ischar(ln), break; end
                ln = strtrim(ln);
                if isempty(ln) || ln(1) == '#' || ln(1) == '%', continue; end
                parts = regexp(ln, '[,;\t ]+', 'split');
                parts = parts(~cellfun(@isempty, parts));
                v = str2double(parts);
                if isempty(v) || any(isnan(v)), continue; end   % header / text row
                rows{end+1,1} = v; %#ok<AGROW>
                nc = max(nc, numel(v));
            end
            A = nan(numel(rows), nc);
            for i = 1:numel(rows), A(i,1:numel(rows{i})) = rows{i}; end
        end

        function [gx, gy, gz] = xyzToGrid(x, y, z)
            ok = isfinite(x) & isfinite(y) & isfinite(z);
            x = x(ok); y = y(ok); z = z(ok);
            ux = unique(x); uy = unique(y);
            if numel(ux)*numel(uy) == numel(x) && numel(ux) > 1 && numel(uy) > 1
                [~, ix] = ismember(x, ux); [~, iy] = ismember(y, uy);
                gz = accumarray([iy ix], z, [numel(uy) numel(ux)], @mean, NaN);
                gx = ux; gy = uy;
            else
                n = min(400, max(50, round(sqrt(numel(x)))));
                gx = linspace(min(x), max(x), n);
                gy = linspace(min(y), max(y), n);
                [XX, YY] = meshgrid(gx, gy);
                try
                    Fs = scatteredInterpolant(x, y, z, 'linear', 'nearest');
                    gz = Fs(XX, YY);
                catch
                    gz = griddata(x, y, z, XX, YY, 'linear');
                end
            end
        end

        function [gx, gy, gz, isGeo] = readNetcdfGrid(file)
            info = ncinfo(file);
            names = {info.Variables.Name};
            vx = OceanEnvironment.pickName(names, {'lon','longitude','x'});
            vy = OceanEnvironment.pickName(names, {'lat','latitude','y'});
            vz = OceanEnvironment.pickName(names, {'elevation','z','depth','bathymetry','topo','band1','height'});
            if isempty(vx) || isempty(vy) || isempty(vz)
                error('Could not find x/y/z variables in %s (found: %s).', file, strjoin(names, ', '));
            end
            gx = double(ncread(file, vx)); gy = double(ncread(file, vy));
            gz = double(ncread(file, vz));
            gx = gx(:); gy = gy(:);
            % orient to ny x nx using the variable's dimension names when possible
            iz = find(strcmp(names, vz), 1);
            dn = {};
            try, dn = lower({info.Variables(iz).Dimensions.Name}); catch, end
            xFirst = numel(dn) >= 2 && any(strcmp(dn{1}, {lower(vx),'lon','longitude','x'}));
            if xFirst || (isempty(dn) && size(gz,1) == numel(gx) && size(gz,2) == numel(gy) && numel(gx) ~= numel(gy))
                gz = gz.';
            end
            isGeo = any(strcmpi(vx, {'lon','longitude'}));
            [gx, gy, gz] = OceanEnvironment.decimateGrid(gx, gy, gz, 1200);
        end

        function [gx, gy, gz, isGeo] = readMatGrid(file)
            Sm = load(file);
            fn = fieldnames(Sm); lf = lower(fn);
            pick = @(c) fn(ismember(lf, c));
            isGeo = [];
            fx = pick({'x','xx','easting'}); fy = pick({'y','yy','northing'});
            flon = pick({'lon','longitude'}); flat = pick({'lat','latitude'});
            fz = pick({'z','depth','elevation','bathy','bathymetry','zz','topo'});
            if ~isempty(flon) && ~isempty(flat), fx = flon; fy = flat; isGeo = true;
            elseif ~isempty(fx), isGeo = false;
            end
            if ~isempty(fx) && ~isempty(fy) && ~isempty(fz)
                X = double(Sm.(fx{1})); Y = double(Sm.(fy{1})); Z = double(Sm.(fz{1}));
                if isvector(X) && isvector(Y) && isequal(size(Z), [numel(Y) numel(X)])
                    gx = X(:); gy = Y(:); gz = Z;
                elseif isvector(X) && isvector(Y) && isequal(size(Z), [numel(X) numel(Y)])
                    gx = X(:); gy = Y(:); gz = Z.';
                elseif isequal(size(X), size(Z)) && isequal(size(Y), size(Z)) && ~isvector(Z)
                    gx = X(1,:); gy = Y(:,1); gz = Z;
                    if abs(gx(end)-gx(1)) < eps   % ndgrid orientation
                        gx = X(:,1); gy = Y(1,:); gz = Z.';
                    end
                else
                    [gx, gy, gz] = OceanEnvironment.xyzToGrid(X(:), Y(:), Z(:));
                end
            else
                A = [];
                for k = 1:numel(fn)
                    v = Sm.(fn{k});
                    if isnumeric(v) && size(v,2) == 3 && size(v,1) > 3, A = double(v); break; end
                end
                if isempty(A), error('MAT file must contain X/Y/Z (or lon/lat/z) or an N x 3 matrix.'); end
                [gx, gy, gz] = OceanEnvironment.xyzToGrid(A(:,1), A(:,2), A(:,3));
            end
            [gx, gy, gz] = OceanEnvironment.decimateGrid(gx, gy, gz, 1200);
        end

        function [gx, gy, gz] = decimateGrid(gx, gy, gz, nmax)
            sx = max(1, ceil(numel(gx)/nmax)); sy = max(1, ceil(numel(gy)/nmax));
            gx = gx(1:sx:end); gy = gy(1:sy:end); gz = gz(1:sy:end, 1:sx:end);
        end

        function nm = pickName(names, candidates)
            nm = '';
            ln = lower(names);
            for k = 1:numel(candidates)
                idx = find(strcmp(ln, candidates{k}), 1);
                if ~isempty(idx), nm = names{idx}; return; end
            end
        end

        function V = readField4(file, info, vname, axisNames, tIdx)
            % Read a lon/lat/depth(/time) variable and return [nlon nlat ndep]
            idx = find(strcmp({info.Variables.Name}, vname), 1);
            dims = info.Variables(idx).Dimensions;
            dn = lower({dims.Name});
            V = double(ncread(file, vname));
            % map dims to lon/lat/dep/time
            order = zeros(1,numel(dn));
            for d = 1:numel(dn)
                if any(strcmpi(dn{d}, {lower(axisNames{1}),'lon','longitude','x'})), order(d) = 1;
                elseif any(strcmpi(dn{d}, {lower(axisNames{2}),'lat','latitude','y'})), order(d) = 2;
                elseif any(strcmpi(dn{d}, {lower(axisNames{3}),'depth','z','lev','level','deptht'})), order(d) = 3;
                else, order(d) = 4;
                end
            end
            tDim = find(order == 4, 1);
            if ~isempty(tDim)
                sz = [size(V) ones(1, max(0, numel(order)-ndims(V)))];   % keep trailing singletons
                tIdx = min(tIdx, sz(tDim));
                subs = repmat({':'}, 1, numel(sz));
                subs{tDim} = tIdx;
                V = V(subs{:});
            end
            keep = find(order ~= 4);
            [~, p] = sort(order(keep));
            V = permute(V, [keep(p) setdiff(1:max(ndims(V), numel(order)), keep(p))]);
            V = reshape(V, size(V,1), size(V,2), size(V,3));
        end

        function y = uniformLookup(zg, vg, z)
            % Fast linear interpolation on a uniform grid zg (column), clamped
            dz = zg(2) - zg(1);
            u = min(max(z(:), zg(1)), zg(end)) / dz;
            i0 = floor(u);
            fr = u - i0;
            i0 = i0 + 1;
            atEnd = i0 >= numel(zg);
            i0(atEnd) = numel(zg) - 1; fr(atEnd) = 1;
            vg = vg(:);
            y = vg(i0) .* (1 - fr) + vg(i0 + 1) .* fr;
            y = reshape(y, size(z));
        end

        function m = nanmedianLocal(v)
            v = v(isfinite(v));
            if isempty(v), m = NaN; else, m = median(v); end
        end
    end
end
