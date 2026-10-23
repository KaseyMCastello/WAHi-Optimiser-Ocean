function cfg = wahi_geometry_cfg(cfg, P, env, Hn, grid)
% wahi_geometry_cfg  Per-geometry localiser settings shared by the GUI and the
% batch runner: search radius, the localiser's sound-speed model and the ray
% model settings (truth propagation and ray-traced localiser).
%
%   cfg  : base localiser config (fs, call, SNR settings, ...)
%   P    : GUI parameter struct (readParams); fields used: locSSPMode, cLoc,
%          propModel, BL, nRays, rayArrivals, dur
%   env  : OceanEnvironment struct (truth)
%   Hn   : N x 3 nominal hydrophone positions (world frame)
%   grid : K x 3 source positions (used for ranges and depths)
%
% The ray model is range-independent: it uses the sound-speed profile at the
% array centroid and a flat bottom at
%   D = max(water depth at the array, deepest hydrophone + 0.5, deepest source + 1).

ctr = mean(Hn(:,1:2), 1);
if isempty(grid), grid = [ctr -10]; end
dmax = max(sqrt(sum((grid(:,1:2) - ctr).^2, 2)));
cfg.searchRmax = 1.25*dmax + 2000;

prof = OceanEnvironment.profileAt(env, ctr(1), ctr(2));
mode = getf(P, 'locSSPMode', 'constant');
switch mode
    case 'env'
        cfg.locSSP = prof; cfg.locModel = 'straight';
    case 'ray'
        cfg.locSSP = prof; cfg.locModel = 'ray';
    otherwise
        cfg.locSSP = OceanEnvironment.constantSSP(getf(P, 'cLoc', 1500)); cfg.locModel = 'straight';
end

wd0 = OceanEnvironment.waterDepth(env, ctr(1), ctr(2));
D = max([wd0, max(-Hn(:,3)) + 0.5, max(-grid(:,3)) + 1, 2]);
D = ceil(D);
rh = 0;
for i = 1:size(Hn,1)
    rh = max(rh, max(sqrt(sum((grid(:,1:2) - Hn(i,1:2)).^2, 2))));
end
o = OceanRayModel.defaults();
o.BL = getf(P, 'BL', o.BL);
o.nRays = getf(P, 'nRays', o.nRays);
o.rmax = 10000*ceil((1.1*rh + 1000)/10000);

cfg.propModel = 'straight';
if strcmp(getf(P, 'propModel', 'straight'), 'ray')
    cfg.propModel = 'ray';
    cfg.ray = struct('ssp', prof, 'D', D, 'o', o, 'arrivals', getf(P, 'rayArrivals', 'all'), ...
        'window', max(1, 3*getf(P, 'dur', 1)));
end
cfg.locRay = struct('ssp', prof, 'D', D, 'o', o, 'nZ', 12);
end

function v = getf(S, name, def)
if isfield(S, name) && ~isempty(S.(name)), v = S.(name); else, v = def; end
end
