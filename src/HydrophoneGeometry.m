% HydrophoneGeometry - array geometry helpers for WAH-i Ocean.
%
% Geometry convention
%   M : N x 3 DESIGN (body-frame) hydrophone coordinates [x y z] in metres,
%       relative to the array origin. Before rotation x = East, y = North,
%       z = up.
%   A : placement struct: ox, oy (array origin, m), mode, refDepth (m),
%       yaw, pitch, roll (deg, optional, default 0)
%
% Orientation (rigid rotation about the array origin, right-hand rule):
%   roll  : about the body y axis  (+roll  tips the +x end DOWN)
%   pitch : about the body x axis  (+pitch tips the +y end UP)
%   yaw   : about the vertical     (+yaw  turns counter-clockwise seen from above,
%                                   i.e. from East towards North)
%   applied in the order roll -> pitch -> yaw:  W = Rz(yaw) Rx(pitch) Ry(roll) M'
%
% Placement modes
%   'seafloor' : every hydrophone sits on the local seafloor; z = height above
%                bottom (HAB). Only yaw is applied (pitch/roll ignored).
%   'frame'    : rigid frame resting on the seafloor at the array origin;
%                world z = -waterDepth(origin) + rotated z.
%   'fixed'    : rigid array at a reference depth; world z = -refDepth + rotated z.
%
% Constraint Lmax is a HORIZONTAL radius about the array centroid.
%
% CSV files written by writeCSV carry '#' metadata lines (origin, mode,
% lat/lon) and a header row; readCSV reads these back, and also accepts
% plain 2- or 3-column numeric files (x y [z]).

classdef HydrophoneGeometry
    methods (Static)

        function names = patternNames()
            names = {'Polygon + centre','Regular polygon','Rectangular grid','Line', ...
                'L-shape','Random (within Lmax)','Square (side dx)', ...
                'Tetrahedron (edge dx)','Cube (edge dx)','Octahedron (edge dx)'};
        end

        function M = pattern(name, o)
            % o: N, dx, dy, nx, ny, rot (deg), jitter (m), z (m), Lmax, minPair
            N = max(1, round(o.N));
            Z = [];                        % body z for 3-D frames (base at 0)
            a = o.dx;
            switch name
                case 'Polygon + centre'
                    k = max(N-1, 1);
                    a = (0:k-1)'*2*pi/k + pi/2;
                    M = [0 0; o.dx*cos(a) o.dx*sin(a)];
                case 'Regular polygon'
                    a = (0:N-1)'*2*pi/N + pi/2;
                    M = [o.dx*cos(a) o.dx*sin(a)];
                case 'Rectangular grid'
                    nx = max(1, round(o.nx)); ny = max(1, round(o.ny));
                    [gx, gy] = meshgrid((0:nx-1)*o.dx, (0:ny-1)*o.dy);
                    M = [gx(:) gy(:)];
                case 'Line'
                    M = [(0:N-1)'*o.dx zeros(N,1)];
                case 'L-shape'
                    nA = ceil((N+1)/2); nB = N - nA + 1;
                    xa = (0:nA-1)'*o.dx;
                    yb = (1:nB-1)'*o.dy;
                    M = [xa zeros(nA,1); zeros(nB-1,1) yb];
                case 'Random (within Lmax)'
                    M = HydrophoneGeometry.random(N, o.Lmax, o.minPair, o.z);
                    M = M(:,1:2);
                case 'Square (side dx)'
                    M = a/2*[-1 -1; 1 -1; 1 1; -1 1];
                case 'Tetrahedron (edge dx)'
                    rb = a/sqrt(3); ang = [90 210 330]'*pi/180;
                    M = [rb*cos(ang) rb*sin(ang); 0 0];
                    Z = [0; 0; 0; a*sqrt(2/3)];
                case 'Cube (edge dx)'
                    [cx, cy, cz] = ndgrid([-1 1]*a/2, [-1 1]*a/2, [0 a]);
                    M = [cx(:) cy(:)]; Z = cz(:);
                case 'Octahedron (edge dx)'
                    hh = a/sqrt(2);
                    M = [hh 0; -hh 0; 0 hh; 0 -hh; 0 0; 0 0];
                    Z = [hh; hh; hh; hh; 0; 2*hh];
                otherwise
                    error('Unknown pattern: %s', name);
            end
            M = M - mean(M,1);
            th = o.rot*pi/180;
            R = [cos(th) -sin(th); sin(th) cos(th)];
            M = (R*M.').';
            if isfield(o,'jitter') && o.jitter > 0
                M = M + o.jitter*randn(size(M));
            end
            if isempty(Z), Z = zeros(size(M,1),1); end
            M = [M Z + o.z];
        end

        function R = orientMatrix(yaw, pitch, roll)
            % Rz(yaw) * Rx(pitch) * Ry(roll), angles in degrees
            y = yaw*pi/180; p = pitch*pi/180; r = roll*pi/180;
            Rz = [cos(y) -sin(y) 0; sin(y) cos(y) 0; 0 0 1];
            Rx = [1 0 0; 0 cos(p) -sin(p); 0 sin(p) cos(p)];
            Ry = [cos(r) 0 sin(r); 0 1 0; -sin(r) 0 cos(r)];
            R = Rz * Rx * Ry;
        end

        function A = fillPlacement(A)
            if ~isfield(A,'yaw') || isempty(A.yaw), A.yaw = 0; end
            if ~isfield(A,'pitch') || isempty(A.pitch), A.pitch = 0; end
            if ~isfield(A,'roll') || isempty(A.roll), A.roll = 0; end
        end

        function W = worldOffsets(M, A)
            % Rotated offsets [East North Up] of each hydrophone from the array
            % origin (seafloor mode: yaw only, z = height above bottom).
            A = HydrophoneGeometry.fillPlacement(A);
            if strcmp(A.mode, 'seafloor')
                R = HydrophoneGeometry.orientMatrix(A.yaw, 0, 0);
            else
                R = HydrophoneGeometry.orientMatrix(A.yaw, A.pitch, A.roll);
            end
            W = (R * M.').';
        end

        function M = random(N, Lmax, minPair, z0)
            if nargin < 4, z0 = 0; end
            for t = 1:5000
                r = Lmax*sqrt(0.0625 + 0.9375*rand(N,1));
                a = 2*pi*rand(N,1);
                M = [r.*cos(a) r.*sin(a) z0*ones(N,1)];
                M(:,1:2) = M(:,1:2) - mean(M(:,1:2),1);
                if HydrophoneGeometry.isValid(M, Lmax, minPair), return; end
            end
            error('Could not generate a valid random geometry (loosen Lmax / minimum spacing).');
        end

        function [ok, msg] = isValid(M, Lmax, minPair)
            ok = true; msg = '';
            if size(M,1) < 3
                ok = false; msg = 'At least 3 hydrophones are required.'; return;
            end
            c = mean(M(:,1:2),1);
            r = sqrt(sum((M(:,1:2) - c).^2, 2));
            if any(r > Lmax*(1+1e-9))
                ok = false; msg = sprintf('A hydrophone is %.0f m from the centroid (limit %.0f m).', max(r), Lmax); return;
            end
            if minPair > 0
                D = HydrophoneGeometry.pairDist(M);
                D(1:size(D,1)+1:end) = inf;
                if any(D(:) < minPair)
                    ok = false; msg = sprintf('Closest pair is %.1f m apart (minimum %.1f m).', min(D(:)), minPair); return;
                end
            end
            Xc = M(:,1:2) - c;
            X3 = M - mean(M,1);
            if rank(X3, 1e-6*max(1,max(abs(X3(:))))) < 2
                ok = false; msg = 'Hydrophones are collinear (no fix possible).'; return;
            end
            if rank(Xc, 1e-6*max(1,max(abs(Xc(:))))) < 2
                msg = 'Note: collinear in plan view (before rotation).';
            end
        end

        function M = recentreXY(M)
            M(:,1:2) = M(:,1:2) - mean(M(:,1:2),1);
        end

        function Mout = shake(M, sigma, Lmax, minPair, moveZ)
            Mout = M;
            Mout(:,1:2) = M(:,1:2) + sigma*randn(size(M,1),2);
            if moveZ, Mout(:,3) = M(:,3) + sigma*randn(size(M,1),1); end
            Mout = HydrophoneGeometry.recentreXY(Mout);
            Mout = HydrophoneGeometry.clampXY(Mout, Lmax);
            if ~HydrophoneGeometry.isValid(Mout, Lmax, minPair), Mout = M; end
        end

        function M = clampXY(M, Lmax)
            r = sqrt(sum(M(:,1:2).^2,2));
            m = r > Lmax;
            if any(m), M(m,1:2) = M(m,1:2) .* (Lmax*0.999 ./ r(m)); end
        end

        function D = pairDist(M)
            n = size(M,1);
            D = zeros(n);
            for i = 1:n
                D(:,i) = sqrt(sum((M - M(i,:)).^2, 2));
            end
        end

        function H = toAbsolute(M, A, env)
            % Absolute hydrophone positions (x, y, z elevation) in the world frame
            W = HydrophoneGeometry.worldOffsets(M, A);
            H = zeros(size(M,1),3);
            H(:,1) = A.ox + W(:,1);
            H(:,2) = A.oy + W(:,2);
            switch A.mode
                case 'seafloor'
                    wd = OceanEnvironment.waterDepth(env, H(:,1), H(:,2));
                    H(:,3) = -wd + max(M(:,3), 0);
                case 'frame'
                    wd0 = OceanEnvironment.waterDepth(env, A.ox, A.oy);
                    H(:,3) = -wd0 + W(:,3);
                otherwise
                    H(:,3) = -A.refDepth + W(:,3);
            end
            H(:,3) = min(H(:,3), -0.5);   % keep below the surface
        end

        function bad = belowSeabed(H, env)
            % Hydrophones deeper than the local seafloor (e.g. a tilted frame)
            wd = OceanEnvironment.waterDepth(env, H(:,1), H(:,2));
            bad = -H(:,3) > wd + 1e-6;
        end

        %% ------------------------------------------------------------------
        function writeCSV(file, M, A, env)
            H = HydrophoneGeometry.toAbsolute(M, A, env);
            hasGeo = isfinite(env.lat0) && isfinite(env.lon0);
            fid = fopen(file, 'w');
            if fid < 0, error('Cannot write %s', file); end
            cleaner = onCleanup(@() fclose(fid));
            fprintf(fid, '# WAH-i Ocean hydrophone geometry\n');
            fprintf(fid, '# origin_x_m=%.3f\n# origin_y_m=%.3f\n', A.ox, A.oy);
            A = HydrophoneGeometry.fillPlacement(A);
            W = HydrophoneGeometry.worldOffsets(M, A);
            fprintf(fid, '# placement_mode=%s\n# reference_depth_m=%.3f\n', A.mode, A.refDepth);
            fprintf(fid, '# yaw_deg=%.4f\n# pitch_deg=%.4f\n# roll_deg=%.4f\n', A.yaw, A.pitch, A.roll);
            if hasGeo
                fprintf(fid, '# origin_lat=%.8f\n# origin_lon=%.8f\n', env.lat0, env.lon0);
            end
            fprintf(fid, '# x/y/z_rel_m = design (body-frame) coordinates; z = height above bottom (seafloor mode) or up from the frame base / reference depth\n');
            fprintf(fid, '# east/north/up_rel_m = offsets from the origin after yaw/pitch/roll\n');
            fprintf(fid, 'id,x_rel_m,y_rel_m,z_rel_m,east_rel_m,north_rel_m,up_rel_m,x_abs_m,y_abs_m,depth_m,lat,lon\n');
            for i = 1:size(M,1)
                if hasGeo
                    [lo, la] = OceanEnvironment.local2geo(H(i,1), H(i,2), env.lon0, env.lat0);
                else
                    lo = NaN; la = NaN;
                end
                fprintf(fid, '%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.3f,%.3f,%.3f,%.8f,%.8f\n', i, M(i,1), M(i,2), M(i,3), ...
                    W(i,1), W(i,2), W(i,3), H(i,1), H(i,2), -H(i,3), la, lo);
            end
        end

        function [M, meta] = readCSV(file)
            % Reads geometry. Returns relative coordinates M (N x 3) and meta
            % struct with any origin/mode metadata found in '#' lines.
            meta = struct('ox',NaN,'oy',NaN,'mode','','refDepth',NaN,'lat0',NaN,'lon0',NaN,'hadHeader',false, ...
                'yaw',NaN,'pitch',NaN,'roll',NaN);
            fid = fopen(file, 'r');
            if fid < 0, error('Cannot open %s', file); end
            cleaner = onCleanup(@() fclose(fid));
            header = {};
            rows = {};
            while true
                ln = fgetl(fid);
                if ~ischar(ln), break; end
                ln = strtrim(ln);
                if isempty(ln), continue; end
                if ln(1) == '#'
                    kv = regexp(ln, '#\s*(\w+)\s*=\s*(\S+)', 'tokens', 'once');
                    if ~isempty(kv)
                        switch lower(kv{1})
                            case 'origin_x_m', meta.ox = str2double(kv{2});
                            case 'origin_y_m', meta.oy = str2double(kv{2});
                            case 'placement_mode', meta.mode = kv{2};
                            case 'reference_depth_m', meta.refDepth = str2double(kv{2});
                            case 'origin_lat', meta.lat0 = str2double(kv{2});
                            case 'origin_lon', meta.lon0 = str2double(kv{2});
                            case 'yaw_deg', meta.yaw = str2double(kv{2});
                            case 'pitch_deg', meta.pitch = str2double(kv{2});
                            case 'roll_deg', meta.roll = str2double(kv{2});
                        end
                    end
                    continue;
                end
                parts = strtrim(regexp(ln, '[,;\t]', 'split'));
                if numel(parts) == 1, parts = regexp(ln, '\s+', 'split'); end
                v = str2double(parts);
                if all(isnan(v)) && isempty(rows)
                    header = lower(parts); meta.hadHeader = true; continue;
                end
                rows{end+1,1} = v; %#ok<AGROW>
            end
            if isempty(rows), error('No numeric rows in %s', file); end
            nc = max(cellfun(@numel, rows));
            A = nan(numel(rows), nc);
            for i = 1:numel(rows), A(i,1:numel(rows{i})) = rows{i}; end
            if ~isempty(header)
                ix = find(strcmp(header,'x_rel_m') | strcmp(header,'x') | strcmp(header,'x_m'), 1);
                iy = find(strcmp(header,'y_rel_m') | strcmp(header,'y') | strcmp(header,'y_m'), 1);
                iz = find(strcmp(header,'z_rel_m') | strcmp(header,'z') | strcmp(header,'z_m'), 1);
                if isempty(ix) || isempty(iy)
                    numCols = find(all(isfinite(A),1));
                    if strcmp(header{1},'id'), numCols = numCols(numCols > 1); end
                    ix = numCols(1); iy = numCols(2);
                    if numel(numCols) >= 3, iz = numCols(3); end
                end
            else
                ix = 1; iy = 2; iz = [];
                if nc >= 3, iz = 3; end
            end
            M = zeros(size(A,1),3);
            M(:,1) = A(:,ix); M(:,2) = A(:,iy);
            if ~isempty(iz), M(:,3) = A(:,iz); end
            M = M(all(isfinite(M),2),:);
            if size(M,1) < 3, error('Need at least 3 hydrophones in %s', file); end
        end

        function q = prctileLocal(v, p)
            v = sort(v(isfinite(v)));
            if isempty(v), q = NaN; return; end
            n = numel(v);
            pos = 1 + (n-1)*p/100;
            lo = floor(pos); hi = ceil(pos);
            q = v(lo) + (pos-lo)*(v(hi)-v(lo));
        end
    end
end
