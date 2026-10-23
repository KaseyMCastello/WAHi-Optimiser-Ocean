function [pts, info] = wahi_build_field_grid(F, env)
% wahi_build_field_grid  Evenly spaced synthetic source positions.
%
%   [pts, info] = wahi_build_field_grid(F, env)
%
% F fields
%   cx, cy     : field centre (m, world frame)
%   R, Rin     : outer / inner horizontal radius (m) about the centre
%   step       : horizontal spacing between sources (m)
%   layout     : 'cartesian' (square grid) or 'polar' (equal-spaced rings,
%                azimuthal spacing ~ step along each ring)
%   azFrom, azTo : compass sector (deg, clockwise from North), 0-360 = all
%   depths     : vector of source depths (m, positive down)
%   seabedClearance : sources must be at least this far above the seabed (m)
%
% Sources on land or below the seabed are dropped.
% pts : K x 3 [x y z] with z = -depth;  info.depthIdx : layer index per point

if ~isfield(F,'seabedClearance'), F.seabedClearance = 1; end
if ~isfield(F,'layout'), F.layout = 'cartesian'; end
depths = abs(F.depths(:))';
if isempty(depths), depths = 10; end
if F.step <= 0, error('Grid spacing must be > 0.'); end
if F.Rin >= F.R, error('Field inner radius must be smaller than the outer radius.'); end

switch lower(F.layout)
    case 'polar'
        xy = [];
        if F.Rin <= 0, xy = [0 0]; r0 = F.step; else, r0 = F.Rin; end
        for r = r0:F.step:F.R + 1e-9
            nA = max(1, round(2*pi*r/F.step));
            a = (0:nA-1)'*2*pi/nA;
            xy = [xy; r*sin(a) r*cos(a)]; %#ok<AGROW>
        end
    otherwise
        n = floor(F.R/F.step);
        ax = (-n:n)*F.step;
        [X, Y] = meshgrid(ax, ax);
        xy = [X(:) Y(:)];
        r = sqrt(sum(xy.^2,2));
        xy = xy(r <= F.R + 1e-9 & r >= F.Rin - 1e-9, :);
end

% compass sector
a1 = mod(F.azFrom, 360); a2 = mod(F.azTo, 360);
if ~(abs(F.azTo - F.azFrom) >= 360 || (a1 == a2 && F.azTo ~= F.azFrom))
    brg = mod(atan2d(xy(:,1), xy(:,2)), 360);
    if a1 <= a2, keep = brg >= a1 & brg <= a2;
    else,        keep = brg >= a1 | brg <= a2;
    end
    keep = keep | (xy(:,1) == 0 & xy(:,2) == 0 & F.Rin <= 0);
    xy = xy(keep,:);
end

xy = xy + [F.cx F.cy];
wd = OceanEnvironment.waterDepth(env, xy(:,1), xy(:,2));

pts = zeros(0,3); depthIdx = zeros(0,1);
nPer = zeros(1, numel(depths));
for k = 1:numel(depths)
    ok = wd - F.seabedClearance >= depths(k) & wd > 0;
    p = [xy(ok,:) -depths(k)*ones(sum(ok),1)];
    pts = [pts; p]; %#ok<AGROW>
    depthIdx = [depthIdx; k*ones(size(p,1),1)]; %#ok<AGROW>
    nPer(k) = size(p,1);
end

info = struct('depths', depths, 'depthIdx', depthIdx, 'nPerDepth', nPer, ...
    'nHorizontal', size(xy,1), 'nDroppedSeabed', numel(depths)*size(xy,1) - size(pts,1));
end
