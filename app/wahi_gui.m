function wahi_gui()
% WAH-i Ocean GUI - Hydrophone Array Design & Optimisation for whale localisation
%
% Adapted from WAH-i (Ravi Umadi, 2026, GPLv3) - originally a microphone-array
% optimiser for echolocating bats. This modified version replaces the in-air
% bat-call model with an underwater model for baleen and toothed whales:
% hydrophones, bathymetry, sound-speed profiles, transmission loss, ambient
% noise, and fields of accuracy out to 100 km.
%
% Tabs
% ----
% 1 Array Design     Build an initial hydrophone layout (patterns with x/y
%                    spacing, or type coordinates into the table), zoomed plan
%                    and side views, pair distances, CSV import/export.
% 2 Environment      Bathymetry (optional), sound speed (constant, 1-D
%                    profile, or 3-D NCOM NetCDF field), sonar-equation SNR,
%                    detection threshold, clock / survey uncertainty.
% 3 Calls            Species presets (fin, blue, humpback, right, minke,
%                    sperm, delphinid, beaked) or an imported template
%                    (waveform samples, frequency contour, or .wav).
% 4 Field & Optimise Evenly spaced synthetic sources in a horizontal field
%                    (up to 100 km), evaluation (WAH), optimisation (WAH-i),
%                    array-position sweep, reports and exports.
%
% Coordinates: x East, y North (m, local), z elevation (m, negative below the
% surface). Hydrophone geometry is stored RELATIVE to the array origin; the
% z column is height above the seafloor (seafloor-mounted) or an offset from
% the reference depth (fixed-depth).
%
% Requirements: MATLAB R2021a or newer (App Designer components). No toolboxes
% are required; the Parallel Computing Toolbox is used if you enable it.
% src/ must be on the path (added automatically when run from the repo).
%
% Licence: GPLv3 (see LICENSE.md). Original WAH-i: https://github.com/raviumadi/WAHi-Optimiser

%% -------------------- Path --------------------
if ~isdeployed
    uHere = fileparts(mfilename('fullpath'));
    uSrc = fullfile(uHere, '..', 'src');
    if exist(uSrc, 'dir') && isempty(strfind(path, uSrc)), addpath(uSrc); end
end
IS_DEPLOYED = isdeployed;

%% -------------------- State --------------------
S = struct();
S.M = [];                       % N x 3 relative hydrophone geometry
S.env = OceanEnvironment.default();
S.sspProfile = [];              % loaded 1-D profile
S.ncom = [];                    % loaded 3-D field
S.callTemplate = struct('values',[],'fs',NaN,'amp',[],'guess','','file','');
S.call = []; S.band = [NaN NaN];
S.gridPts = []; S.gridInfo = [];
S.stats = []; S.diag = struct();
S.H = []; S.runInfo = [];
S.sweep = [];
S.stopFlag = false;
S.busy = false;
S.selRows = [];
S.progressPrefix = '';
S.lastEnvDisplay = '';

CTRL = struct();    % key -> UI handle
PARAMCB = @(src,~) onParamChanged(src.Tag);   % shared ValueChangedFcn (key stored in Tag)
PROGRESSCB = @(i,n) progressCb(i,n);
DEF = struct();     % key -> default value
GRP = struct();     % key -> group name

THEME = struct();
THEME.baseBg       = [0.97 0.985 0.995];
THEME.panelBg      = [0.995 1.0 1.0];
THEME.mutedBg      = [0.97 0.98 1.00];
THEME.panelAlt     = [1.00 0.995 0.98];
THEME.textColor    = [0.10 0.12 0.15];
THEME.accentColor  = [0.00 0.45 0.70];
THEME.secondaryBtn = [0.88 0.93 0.97];
THEME.headerTeal   = [0.62 0.80 0.86];

%% -------------------- Figure --------------------
uScr = get(0, 'ScreenSize');
uW = max(1200, min(1600, uScr(3) - 60)); uH = max(700, min(900, uScr(4) - 90));
fig = uifigure('Name','WAH-i Ocean - Hydrophone Array Design & Optimisation', ...
    'Position',[30 40 uW uH]);
fig.Color = THEME.baseBg;

uMenu = uimenu(fig, 'Text', 'Menu');
uimenu(uMenu, 'Text', 'Save session...', 'MenuSelectedFcn', @(~,~)doSaveSession());
uimenu(uMenu, 'Text', 'Load session...', 'MenuSelectedFcn', @(~,~)doLoadSession());
uimenu(uMenu, 'Text', 'Reset defaults', 'MenuSelectedFcn', @(~,~)doResetDefaults(), 'Separator','on');
uimenu(uMenu, 'Text', 'Tips', 'MenuSelectedFcn', @(~,~)showTips(), 'Separator','on');
uimenu(uMenu, 'Text', 'About', 'MenuSelectedFcn', @(~,~)showAbout());

uRoot = uigridlayout(fig, [2 1]);
uRoot.RowHeight = {'1x', 24};
uRoot.Padding = [8 6 8 6];
uRoot.RowSpacing = 4;

tg = uitabgroup(uRoot);
tg.Layout.Row = 1;
tabArr  = uitab(tg, 'Title', '1  Array Design');
tabEnv  = uitab(tg, 'Title', '2  Environment');
tabCall = uitab(tg, 'Title', '3  Calls');
tabFoA  = uitab(tg, 'Title', '4  Field & Optimise');

lblStatus = uilabel(uRoot, 'Text', 'Ready', 'FontWeight', 'bold');
lblStatus.Layout.Row = 2;
lblStatus.FontColor = THEME.textColor;

%% ====================================================================
%  TAB 1 - ARRAY DESIGN
%% ====================================================================
gArr = uigridlayout(tabArr, [1 2]);
gArr.ColumnWidth = {470, '1x'};
gArr.Padding = [8 8 8 8];
gArr.ColumnSpacing = 10;

gArrL = uigridlayout(gArr, [3 1]);
gArrL.Layout.Row = 1; gArrL.Layout.Column = 1;
gArrL.RowHeight = {300, 292, '1x'};
gArrL.Scrollable = 'on';
gArrL.Padding = [0 0 0 0];
gArrL.RowSpacing = 8;

gPlace = mkPanelGrid(gArrL, 'Placement, orientation & constraints', 11, 1, 1, THEME.panelBg);
addNum(gPlace, 'Array origin East X (m):', 'ox', 0, [-2e5 2e5], '%.1f', ...
    'Array origin in the local frame (m East of the map origin). Hydrophone coordinates are relative to this point.', false, 'array');
addNum(gPlace, 'Array origin North Y (m):', 'oy', 0, [-2e5 2e5], '%.1f', ...
    'Array origin in the local frame (m North of the map origin).', false, 'array');
addDrop(gPlace, 'Placement:', 'placeMode', {'Seafloor-mounted (each element)','Seafloor frame (rigid)','Fixed depth (rigid)'}, ...
    'Seafloor-mounted (each element)', ['Each element: every hydrophone sits on the local seafloor, z = height above bottom (yaw only). ' ...
    'Seafloor frame: rigid frame whose origin rests on the seafloor; z = height above the frame base; yaw/pitch/roll applied. ' ...
    'Fixed depth: rigid array at the reference depth; z = up from the reference point; yaw/pitch/roll applied.'], 'array');
addNum(gPlace, 'Reference depth (m) [fixed mode]:', 'refDepth', 500, [0 11000], '%.1f', ...
    'Depth of the array origin when Placement = Fixed depth.', false, 'array');
addNum(gPlace, 'Yaw (deg, CCW from East):', 'yaw', 0, [-360 360], '%.2f', ...
    'Rotation about the vertical through the array origin. +90 turns the body +x axis from East to North.', false, 'array');
addNum(gPlace, 'Pitch (deg, about body x):', 'pitch', 0, [-180 180], '%.2f', ...
    'Tilt about the body x axis; + raises the +y end. Applied before yaw. Ignored for each-element seafloor placement.', false, 'array');
addNum(gPlace, 'Roll (deg, about body y):', 'roll', 0, [-180 180], '%.2f', ...
    'Tilt about the body y axis; + lowers the +x end. Applied first (roll -> pitch -> yaw). Ignored for each-element seafloor placement.', false, 'array');
addNum(gPlace, 'Max array radius Lmax (m):', 'Lmax', 2000, [1 1e5], '%.1f', ...
    'Maximum HORIZONTAL distance of any hydrophone from the array centroid (optimiser constraint).', false, 'array');
addNum(gPlace, 'Min hydrophone spacing (m):', 'minPair', 100, [0 1e5], '%.1f', ...
    'Minimum allowed distance between any pair of hydrophones (optimiser constraint).', false, 'array');
lblArrDepth = addInfo(gPlace, 'Seafloor at origin: ---', 2);

gPat = mkPanelGrid(gArrL, 'Pattern generator', 10, 2, 1, THEME.mutedBg);
addDrop(gPat, 'Pattern:', 'patName', HydrophoneGeometry.patternNames(), 'Polygon + centre', ...
    'Starting layout. Use the table afterwards for uneven / custom spacing.', 'none');
addNum(gPat, 'Number of hydrophones:', 'patN', 4, [3 64], '%d', 'Used by polygon, line, L-shape and random patterns.', true, 'none');
addNum(gPat, 'Spacing dx (m):', 'patDx', 1000, [0.01 1e5], '%g', ...
    'Polygon radius (centre-to-vertex), grid/line x spacing, L-shape x-arm spacing, or the side/edge length of the square, tetrahedron, cube and octahedron.', false, 'none');
addNum(gPat, 'Spacing dy (m):', 'patDy', 1000, [0.01 1e5], '%g', 'Grid y spacing, or L-shape y-arm spacing.', false, 'none');
addNum(gPat, 'Grid columns nx:', 'patNx', 2, [1 32], '%d', 'Rectangular grid only.', true, 'none');
addNum(gPat, 'Grid rows ny:', 'patNy', 2, [1 32], '%d', 'Rectangular grid only.', true, 'none');
addNum(gPat, 'Rotation (deg, CCW):', 'patRot', 0, [-360 360], '%.1f', 'Rotates the pattern counter-clockwise about its centroid.', false, 'none');
addNum(gPat, 'Random jitter sigma (m):', 'patJit', 0, [0 1e4], '%.1f', 'Adds Gaussian jitter to x/y so the spacing is not perfectly even.', false, 'none');
addNum(gPat, 'z: height above bottom / offset (m):', 'patZ', 5, [-5000 5000], '%.1f', ...
    'Seafloor mode: height above bottom. Fixed mode: vertical offset from reference depth (positive = shallower).', false, 'none');
addButtons(gPat, {'Generate pattern'}, {@(~,~)doGeneratePattern()}, true);

pTbl = uipanel(gArrL, 'Title', 'Hydrophones (click a cell to edit)', 'FontWeight', 'bold');
pTbl.Layout.Row = 3; pTbl.BackgroundColor = THEME.panelBg;
gTbl = uigridlayout(pTbl, [3 1]);
gTbl.RowHeight = {'1x', 28, 28};
gTbl.Padding = [6 6 6 6];
gTbl.RowSpacing = 6;
tblHyd = uitable(gTbl, 'Data', zeros(0,4), 'ColumnName', {'#','x (m)','y (m)','z HAB (m)'}, ...
    'ColumnEditable', [false true true true], 'RowName', {}, ...
    'CellEditCallback', @(src,evt)onTableEdit(src,evt), ...
    'CellSelectionCallback', @(src,evt)onTableSelect(src,evt));
tblHyd.Layout.Row = 1;
uB = addButtons(gTbl, {'Add','Remove selected','Re-centre','Sort by #'}, ...
    {@(~,~)doAddHyd(), @(~,~)doRemoveHyd(), @(~,~)doRecentre(), @(~,~)doRenumber()}, false);
uB.Layout.Row = 2;
uB = addButtons(gTbl, {'Import CSV...','Export CSV...','Evaluate this array (tab 4)'}, ...
    {@(~,~)doImportGeometry(), @(~,~)doExportGeometry(), @(~,~)goEvaluate()}, false);
uB.Layout.Row = 3;

gArrR = uigridlayout(gArr, [3 2]);
gArrR.Layout.Row = 1; gArrR.Layout.Column = 2;
gArrR.RowHeight = {26, '2x', '1x'};
gArrR.ColumnWidth = {'1.3x', '1x'};
gArrR.Padding = [0 0 0 0];
gArrR.RowSpacing = 6; gArrR.ColumnSpacing = 8;
gArrTools = uigridlayout(gArrR, [1 4]);
gArrTools.Layout.Row = 1; gArrTools.Layout.Column = [1 2];
gArrTools.Padding = [0 0 0 0];
gArrTools.ColumnWidth = {'fit','fit','fit','1x'};
cbPairLines = uicheckbox(gArrTools, 'Text', 'Show pair distances', 'Value', true, ...
    'ValueChangedFcn', @(~,~)refreshArrayTab());
cbArrBathy = uicheckbox(gArrTools, 'Text', 'Show bathymetry contours', 'Value', true, ...
    'ValueChangedFcn', @(~,~)refreshArrayTab());
cbArrLmax = uicheckbox(gArrTools, 'Text', 'Show Lmax circle', 'Value', true, ...
    'ValueChangedFcn', @(~,~)refreshArrayTab());
lblArrValid = uilabel(gArrTools, 'Text', '', 'FontWeight', 'bold');

axArrPlan = uiaxes(gArrR);
axArrPlan.Layout.Row = 2; axArrPlan.Layout.Column = [1 2];
axArrSide = uiaxes(gArrR);
axArrSide.Layout.Row = 3; axArrSide.Layout.Column = 1;
gArrInfo = uigridlayout(gArrR, [2 1]);
gArrInfo.Layout.Row = 3; gArrInfo.Layout.Column = 2;
gArrInfo.RowHeight = {'fit', '1x'};
gArrInfo.Padding = [0 0 0 0];
lblArrStats = uilabel(gArrInfo, 'Text', '', 'WordWrap', 'on', 'VerticalAlignment', 'top');
uTgInfo = uitabgroup(gArrInfo);
uTabW = uitab(uTgInfo, 'Title', 'World coordinates');
uTabP = uitab(uTgInfo, 'Title', 'Pair distances');
uGW = uigridlayout(uTabW, [1 1]); uGW.Padding = [2 2 2 2];
tblWorld = uitable(uGW, 'Data', [], 'RowName', {}, ...
    'ColumnName', {'#','East (m)','North (m)','Up (m)','Depth (m)'});
uGP = uigridlayout(uTabP, [1 1]); uGP.Padding = [2 2 2 2];
tblPair = uitable(uGP, 'Data', [], 'RowName', {});

%% ====================================================================
%  TAB 2 - ENVIRONMENT
%% ====================================================================
gEnv = uigridlayout(tabEnv, [1 2]);
gEnv.ColumnWidth = {470, '1x'};
gEnv.Padding = [8 8 8 8];
gEnv.ColumnSpacing = 10;

gEnvL = uigridlayout(gEnv, [4 1]);
gEnvL.Layout.Row = 1; gEnvL.Layout.Column = 1;
gEnvL.RowHeight = {216, 268, 372, 112};
gEnvL.Padding = [0 0 0 0];
gEnvL.RowSpacing = 8;
gEnvL.Scrollable = 'on';

gBathy = mkPanelGrid(gEnvL, 'Bathymetry (optional)', 7, 1, 1, THEME.panelBg);
lblBathy = addInfo(gBathy, 'No bathymetry loaded - flat bottom.', 2);
addButtons(gBathy, {'Load bathymetry...','Flip depth sign','Clear'}, ...
    {@(~,~)doLoadBathy(), @(~,~)doFlipBathy(), @(~,~)doClearBathy()}, true);
addNum(gBathy, 'Flat-bottom water depth (m):', 'flatDepth', 1000, [1 11000], '%.1f', ...
    'Water depth used everywhere when no bathymetry file is loaded.', false, 'env');
addTxt(gBathy, 'Origin latitude (deg):', 'lat0', '', ...
    'Geographic origin of the local x/y frame. Leave blank to use the centre of the first geographic file loaded.', 'env');
addTxt(gBathy, 'Origin longitude (deg):', 'lon0', '', ...
    'Geographic origin of the local x/y frame (decimal degrees, East positive).', 'env');

gSSP = mkPanelGrid(gEnvL, 'Sound speed', 9, 2, 1, THEME.mutedBg);
addDrop(gSSP, 'Truth model (simulation):', 'sspModel', {'Constant','1-D profile (file)','3-D field (NCOM NetCDF)'}, 'Constant', ...
    'Sound speed used to SIMULATE arrivals (the real ocean).', 'env');
addNum(gSSP, 'Constant c (m/s):', 'c0', 1500, [1300 1700], '%.1f', 'Used when the truth model is Constant.', false, 'env');
addButtons(gSSP, {'Load 1-D profile...','Load NCOM file...'}, {@(~,~)doLoadProfile(), @(~,~)doLoadNCOM()}, true);
addNum(gSSP, 'NCOM time index:', 'ncomT', 1, [1 1e5], '%d', 'Time slice used from the NetCDF file (re-load after changing).', true, 'none');
lblSSP = addInfo(gSSP, 'Constant sound speed.', 2);
addDrop(gSSP, 'Localiser assumes:', 'locSSPMode', {'Constant c (below)','Environment profile at array','Ray-traced (profile at array)'}, 'Constant c (below)', ...
    ['Sound speed / propagation the LOCALISER uses. A mismatch with the truth model shows up as localisation error. ' ...
    'Ray-traced: the solver uses first-arrival times from rays traced through the profile at the array (accounts for refraction).'], 'loc');
addNum(gSSP, 'Localiser c (m/s):', 'cLoc', 1500, [1300 1700], '%.1f', 'Constant sound speed assumed by the localiser.', false, 'loc');

gProp = mkPanelGrid(gEnvL, 'Propagation, noise & detection', 13, 3, 1, THEME.panelBg);
addDrop(gProp, 'Propagation model (truth):', 'propModel', {'Straight rays (fast)','Ray tracing (refraction, ducting, multipath)'}, ...
    'Straight rays (fast)', ['Straight rays: travel time along the straight path through the sound-speed profile; TL from the spreading law. ' ...
    'Ray tracing: Snell''s-law rays through the profile at the array with surface/bottom reflections; travel times, TL, shadow zones ' ...
    'and multipath come from the eigenrays.'], 'env');
addNum(gProp, 'Bottom loss per bounce (dB) [rays]:', 'BL', 6, [0 100], '%.1f', ...
    'Reflection loss applied each time a ray reflects off the seafloor (depends on sediment and grazing angle; 0-15 dB typical at low frequency).', false, 'env');
addDrop(gProp, 'Received arrivals [rays]:', 'rayArrivals', {'All eigenrays (multipath)','Strongest arrival only','First arrival only'}, ...
    'All eigenrays (multipath)', ['All: every eigenray within 30 dB of the strongest is summed into the signal; SNR from the total received level. ' ...
    'First / Strongest: the signal contains only that arrival and its SNR comes from that arrival alone (other paths treated as absent). ' ...
    'Multipath can make the cross-correlation lock onto different paths on different hydrophones.'], 'env');
addNum(gProp, 'Rays per launch direction:', 'nRays', 400, [50 5000], '%d', ...
    'Angular resolution of the ray fan (0.05-89 deg, up and down). More rays = finer eigenrays, slower.', true, 'env');
addDrop(gProp, 'SNR model:', 'snrMode', {'Sonar equation','Fixed SNR'}, 'Sonar equation', ...
    'Sonar equation: SNR = SL - TL(r,f) - NL per hydrophone. Fixed: same SNR at every hydrophone and range.', 'env');
addNum(gProp, 'Fixed SNR (dB):', 'snrFixed', 20, [-30 100], '%.1f', 'Used when SNR model = Fixed SNR.', false, 'env');
addNum(gProp, 'Noise level NL (dB re 1 uPa, in band):', 'NL', 88, [0 200], '%.1f', ...
    'Ambient noise level integrated over the call band (spectrum level + 10log10(bandwidth)).', false, 'env');
addDrop(gProp, 'Spreading [straight rays]:', 'spreading', {'Hybrid (spherical then cylindrical)','Spherical (20 log r)','Practical (15 log r)','Cylindrical (10 log r)'}, ...
    'Hybrid (spherical then cylindrical)', 'Geometric spreading law for transmission loss.', 'env');
addNum(gProp, 'Transition range (m) [hybrid]:', 'transR', 1000, [1 1e5], '%.0f', ...
    'Range where spreading changes from spherical to cylindrical (typically ~ water depth).', false, 'env');
addChk(gProp, 'Include seawater absorption (Thorp)', 'absorption', true, 'Frequency-dependent absorption at the call centre frequency.', 'env');
addNum(gProp, 'Detection threshold DT (dB):', 'DT', 10, [-30 100], '%.1f', ...
    'A hydrophone only contributes a TDOA if its received SNR >= DT.', false, 'env');
addNum(gProp, 'Min hydrophones detecting:', 'minDetect', 4, [3 64], '%d', ...
    'Sources heard on fewer hydrophones are "undetected" (2-D solve needs >= 3, 3-D needs >= 4).', true, 'env');

gUnc = mkPanelGrid(gEnvL, 'Uncertainties', 3, 4, 1, THEME.mutedBg);
addNum(gUnc, 'Clock / timing error sigma (ms):', 'timingMs', 0, [0 1e4], '%.3f', ...
    'Random timing error added to every arrival (clock drift / sync error between recorders).', false, 'env');
addNum(gUnc, 'Hydrophone position error sigma (m):', 'posSigma', 0, [0 1e4], '%.2f', ...
    'Survey error: true positions differ from the positions the localiser uses.', false, 'env');

gEnvR = uigridlayout(gEnv, [3 2]);
gEnvR.Layout.Row = 1; gEnvR.Layout.Column = 2;
gEnvR.RowHeight = {'1.3x', '1x', '1x'};
gEnvR.ColumnWidth = {'2x', '1x'};
gEnvR.Padding = [0 0 0 0];
axBathy = uiaxes(gEnvR); axBathy.Layout.Row = 1; axBathy.Layout.Column = 1;
axSSP = uiaxes(gEnvR);   axSSP.Layout.Row = 1;   axSSP.Layout.Column = 2;
axSNR = uiaxes(gEnvR);   axSNR.Layout.Row = 2;   axSNR.Layout.Column = [1 2];
axRays = uiaxes(gEnvR);  axRays.Layout.Row = 3;  axRays.Layout.Column = [1 2];

%% ====================================================================
%  TAB 3 - CALLS
%% ====================================================================
gCall = uigridlayout(tabCall, [1 2]);
gCall.ColumnWidth = {470, '1x'};
gCall.Padding = [8 8 8 8];
gCall.ColumnSpacing = 10;

gCallL = uigridlayout(gCall, [2 1]);
gCallL.Layout.Row = 1; gCallL.Layout.Column = 1;
gCallL.RowHeight = {432, '1x'};
gCallL.Scrollable = 'on';
gCallL.Padding = [0 0 0 0];

gCP = mkPanelGrid(gCallL, 'Call synthesis', 15, 1, 1, THEME.panelBg);
addDrop(gCP, 'Preset:', 'callPreset', WhaleCallLibrary.presetNames(), 'Fin whale - 20 Hz pulse', ...
    'Species preset (approximate). Choosing a preset fills in the fields below; you can then edit them.', 'preset');
lblCallNote = addInfo(gCP, '', 2);
addNum(gCP, 'Sample rate fs (Hz):', 'fs', 500, [20 1e6], '%.0f', ...
    'Simulation sample rate. Must exceed twice the highest call frequency. Low rates keep long-range runs fast.', false, 'call');
addNum(gCP, 'Call duration (s):', 'dur', 1, [1e-5 300], '%.4g', 'Length of the synthetic call (also used to stretch a frequency-contour template).', false, 'call');
addNum(gCP, 'Start frequency (Hz):', 'f0', 23, [0.1 5e5], '%.6g', 'FM/tonal: start frequency. Pulse: lower band edge.', false, 'call');
addNum(gCP, 'End frequency (Hz):', 'f1', 17, [0.1 5e5], '%.6g', 'FM: end frequency. Pulse: upper band edge. Ignored for tonal.', false, 'call');
addDrop(gCP, 'Sweep shape:', 'sweep', {'linear','exponential'}, 'linear', 'Shape of the frequency sweep for FM calls.', 'call');
addNum(gCP, 'Taper (%):', 'taper', 10, [0 100], '%.0f', 'Tukey window taper applied to the call envelope.', false, 'call');
addNum(gCP, 'Source level SL (dB re 1 uPa @ 1 m):', 'SL', 189, [60 260], '%.1f', ...
    'Used by the sonar equation (Environment tab) to compute received SNR.', false, 'env');
addNum(gCP, 'Template sample rate (Hz, 0 = fs):', 'tplFs', 0, [0 1e6], '%.0f', ...
    'Sample rate of an imported WAVEFORM template (ignored for .wav files, which carry their own rate).', false, 'call');
addButtons(gCP, {'Import template...','Refresh preview'}, {@(~,~)doImportTemplate(), @(~,~)updateCall(true)}, true);
lblTpl = addInfo(gCP, 'No template loaded.', 2);
lblBand = addInfo(gCP, 'Band: ---', 1);

pCallHelp = uipanel(gCallL, 'Title', 'Template formats', 'FontWeight', 'bold');
pCallHelp.Layout.Row = 2; pCallHelp.BackgroundColor = THEME.mutedBg;
gCH = uigridlayout(pCallHelp, [1 1]);
uilabel(gCH, 'WordWrap', 'on', 'VerticalAlignment', 'top', 'Text', ...
    ['Import a single vector (e.g. 200 points) from .csv/.txt/.mat:' newline ...
    '  - Frequency contour: values in Hz, spread evenly over the Call duration.' newline ...
    '    Optional 2nd column = relative amplitude.' newline ...
    '  - Waveform: pressure samples at the Template sample rate.' newline ...
    'Or import a .wav clip of a real call (its own sample rate is used).' newline ...
    'You will be asked which interpretation to use.']);

gCallR = uigridlayout(gCall, [2 1]);
gCallR.Layout.Row = 1; gCallR.Layout.Column = 2;
gCallR.Padding = [0 0 0 0];
axWave = uiaxes(gCallR); axWave.Layout.Row = 1;
axSpec = uiaxes(gCallR); axSpec.Layout.Row = 2;

%% ====================================================================
%  TAB 4 - FIELD & OPTIMISE
%% ====================================================================
gFoA = uigridlayout(tabFoA, [1 2]);
gFoA.ColumnWidth = {800, '1x'};
gFoA.Padding = [8 8 8 8];
gFoA.ColumnSpacing = 10;

gFL = uigridlayout(gFoA, [3 2]);
gFL.Layout.Row = 1; gFL.Layout.Column = 1;
gFL.RowHeight = {346, 346, 270};
gFL.Scrollable = 'on';
gFL.ColumnWidth = {'1x', '1x'};
gFL.Padding = [0 0 0 0];
gFL.RowSpacing = 8; gFL.ColumnSpacing = 8;

gFld = mkPanelGrid(gFL, 'Field of accuracy (synthetic sources)', 12, 1, 1, THEME.panelBg);
addNum(gFld, 'Field centre East X (m):', 'fcx', 0, [-2e5 2e5], '%.1f', 'Centre of the study area (fixed in the world frame).', false, 'field');
addNum(gFld, 'Field centre North Y (m):', 'fcy', 0, [-2e5 2e5], '%.1f', 'Centre of the study area.', false, 'field');
addNum(gFld, 'Outer radius R (m):', 'R', 20000, [10 1e5], '%.0f', 'Sources are placed out to this horizontal distance from the field centre (max 100 km).', false, 'field');
addNum(gFld, 'Inner radius Rin (m):', 'Rin', 0, [0 1e5], '%.0f', 'No sources inside this radius (0 = fill the disc).', false, 'field');
addNum(gFld, 'Source spacing (m):', 'step', 2000, [1 5e4], '%.0f', ...
    'Horizontal distance between neighbouring synthetic sources. Smaller = denser and slower.', false, 'field');
addDrop(gFld, 'Layout:', 'layout', {'Cartesian grid','Polar rings'}, 'Cartesian grid', ...
    'Cartesian: square grid. Polar: rings every "spacing" with ~"spacing" between sources along each ring.', 'field');
addNum(gFld, 'Sector from (deg true):', 'azFrom', 0, [-360 720], '%.1f', 'Compass sector (clockwise from North). 0 to 360 = all directions.', false, 'field');
addNum(gFld, 'Sector to (deg true):', 'azTo', 360, [-360 720], '%.1f', 'End of the compass sector (clockwise).', false, 'field');
addTxt(gFld, 'Source depth(s) (m):', 'depths', '50', 'One depth or a comma-separated list, e.g. "20, 100, 500".', 'field');
addNum(gFld, 'Min source-hydrophone dist (m):', 'minDistAny', 100, [0 1e5], '%.1f', 'Sources closer than this to any hydrophone are skipped.', false, 'field');
lblGrid = addInfo(gFld, 'Sources: ---', 2);

gLoc = mkPanelGrid(gFL, 'Localisation & scoring', 12, 2, 1, THEME.mutedBg);
addDrop(gLoc, 'Solve for:', 'solveMode', {'2-D (x,y) at assumed depth','3-D (x,y,z)'}, '2-D (x,y) at assumed depth', ...
    '2-D is typical for near-planar seafloor arrays; 3-D needs depth aperture and >= 4 detecting hydrophones.', 'loc');
addNum(gLoc, 'Assumed source depth (m) [2-D]:', 'assumedDepth', 50, [0 11000], '%.1f', ...
    'Depth assumed by the 2-D solver when "known depth" is off. Set it different from the true depth to test depth-mismatch error.', false, 'loc');
addChk(gLoc, '2-D solver knows each source''s true depth', 'useTrueDepth', true, ...
    'On: the 2-D solve uses the true source depth (isolates horizontal geometry). Off: uses the assumed depth above for every source.', 'loc');
addDrop(gLoc, 'TDOA method:', 'tdoaMethod', {'Waveform cross-correlation','Envelope cross-correlation'}, 'Waveform cross-correlation', ...
    'Envelope correlation is more robust for narrowband tonal calls but less precise.', 'loc');
addDrop(gLoc, 'Score on:', 'scoreMetric', {'Horizontal error','3-D error'}, 'Horizontal error', 'Error used for pass/fail and the optimiser score.', 'loc');
addDrop(gLoc, 'Undetected sources:', 'undetPolicy', {'Count as failures','Exclude from score'}, 'Count as failures', ...
    'Whether sources heard on too few hydrophones count against the pass rate.', 'loc');
addNum(gLoc, 'Target error (m):', 'best', 100, [0.01 1e6], '%.2f', 'Desired localisation error.', false, 'loc');
addNum(gLoc, 'Tolerance margin (m):', 'tol', 50, [0 1e6], '%.2f', 'Pass threshold = Target error + Tolerance.', false, 'loc');
addNum(gLoc, 'Mean-error weight (lambda):', 'lambda', 0.15, [0 10], '%.3f', 'J = passRate - lambda * min(meanErr / penalty, 1).', false, 'loc');
addNum(gLoc, 'Failure penalty (m):', 'failPen', 5000, [1 1e7], '%.0f', 'Error assigned to failed/undetected sources (also caps errors).', false, 'loc');
addNum(gLoc, 'Random seed:', 'seed', 1, [0 2^31-2], '%d', 'Seeds noise, clock and position-error realisations (same noise for every candidate geometry).', true, 'loc');
addChk(gLoc, 'Use parallel pool (Parallel Computing Toolbox)', 'useParallel', false, 'Evaluates sources in parallel with parfor. No progress updates while running.', 'none');

gOpt = mkPanelGrid(gFL, 'Optimiser (WAH-i)', 9, 3, 1, THEME.panelBg);
addNum(gOpt, 'Max iterations:', 'maxIters', 25, [1 5000], '%d', 'One hydrophone is perturbed per iteration.', true, 'opt');
addNum(gOpt, 'Probes per hydrophone:', 'probes', 8, [1 16], '%d', 'Number of directions tried per iteration (plus 2 radial).', true, 'opt');
addNum(gOpt, 'Initial step fraction:', 'initStep', 0.25, [0 1], '%.3f', 'Step = fraction x Lmax, decreasing linearly over iterations.', false, 'opt');
addNum(gOpt, 'Minimum step fraction:', 'minStep', 0.03, [0 1], '%.3f', 'Lower bound on the step fraction.', false, 'opt');
addNum(gOpt, 'Stall limit:', 'stall', 12, [1 1e4], '%d', 'Non-improving iterations before a shake/reset.', true, 'opt');
addNum(gOpt, 'Shake fraction:', 'shake', 0.08, [0 1], '%.3f', 'Shake magnitude as a fraction of Lmax.', false, 'opt');
addNum(gOpt, 'Min improvement (dJ):', 'epsJ', 1e-6, [0 1], '%.1e', 'Minimum score gain to accept a move.', false, 'opt');
addNum(gOpt, 'Target pass rate:', 'targetPass', 0.95, [0 1], '%.2f', 'Stop early at this pass rate.', false, 'opt');
addDrop(gOpt, 'Moves:', 'moveMode', {'Horizontal (x,y) only','3-D (x,y,z)'}, 'Horizontal (x,y) only', ...
    'Seafloor arrays normally only move horizontally.', 'opt');

pAct = uipanel(gFL, 'Title', 'Workflow', 'FontWeight', 'bold');
pAct.Layout.Row = 1; pAct.Layout.Column = 2; pAct.BackgroundColor = THEME.panelAlt;
gAct = uigridlayout(pAct, [6 3]);
gAct.RowHeight = {34, 34, 30, 30, 30, '1x'};
gAct.ColumnWidth = {'1x','1x','1x'};
gAct.Padding = [8 8 8 8]; gAct.RowSpacing = 6; gAct.ColumnSpacing = 6;
btnRun = mkBtn(gAct, '1  Evaluate field (WAH)', @(~,~)doRunWAH(), 1, [1 3], true, ...
    'Simulate every synthetic source with the current array and environment.');
btnOpt = mkBtn(gAct, '2  Optimise geometry (WAH-i)', @(~,~)doOptimiser(), 2, [1 3], true, ...
    'Iteratively move hydrophones to improve the pass rate.');
btnStop = mkBtn(gAct, 'Stop', @(~,~)doStop(), 3, 1, false, 'Stop the running evaluation / optimisation / sweep.');
btnStop.Enable = 'off';
btnReport = mkBtn(gAct, 'Report', @(~,~)doReport(), 3, 2, false, 'Before/after report for the last optimisation run.');
btnReport.Enable = 'off';
btnExportRun = mkBtn(gAct, 'Export run (.mat)', @(~,~)doExportRun(), 3, 3, false, 'Save the run history H and runInfo.');
btnExportRun.Enable = 'off';
mkBtn(gAct, 'Export results CSV', @(~,~)doExportResults(), 4, 1, false, 'Per-source errors, SNR and detections of the last evaluation.');
mkBtn(gAct, 'Inspect a source...', @(~,~)doInspect(), 4, 2, false, 'Simulate one source and view the received signals and estimate.');
btnExportWS = mkBtn(gAct, 'Export workspace', @(~,~)doExportWorkspace(), 4, 3, false, 'Assign results to the MATLAB base workspace.');
if IS_DEPLOYED, btnExportWS.Enable = 'off'; end
mkBtn(gAct, 'Save session...', @(~,~)doSaveSession(), 5, 1, false, 'Save every setting, the array, bathymetry, SSP and template.');
mkBtn(gAct, 'Load session...', @(~,~)doLoadSession(), 5, 2, false, 'Restore a saved session.');
mkBtn(gAct, 'Export geometry CSV...', @(~,~)doExportGeometry(), 5, 3, false, 'Same as Export CSV on the Array Design tab.');
lblStats = uilabel(gAct, 'Text', 'No results yet.', 'WordWrap', 'on', 'VerticalAlignment', 'top', 'FontSize', 12);
lblStats.Layout.Row = 6; lblStats.Layout.Column = [1 3];

gSw = mkPanelGrid(gFL, 'Array position & orientation sweep', 8, 2, 2, THEME.mutedBg);
addDrop(gSw, 'Candidate origins:', 'sweepArea', {'Within the field (radius R)','Around current origin'}, 'Within the field (radius R)', ...
    'Where to try placing the array (same geometry, translated).', 'none');
addNum(gSw, 'Sweep radius (m) [around origin]:', 'sweepR', 10000, [0 1e5], '%.0f', 'Radius of candidate origins around the current origin.', false, 'none');
addNum(gSw, 'Sweep spacing (m):', 'sweepStep', 5000, [1 1e5], '%.0f', 'Spacing between candidate array origins.', false, 'none');
addTxt(gSw, 'Yaw angles to try (deg):', 'sweepYaw', '', ...
    'Optional comma-separated yaw angles, e.g. "0, 30, 60". Blank = keep the current yaw. Each position is evaluated at every yaw.', 'none');
addDrop(gSw, 'Rank by:', 'sweepRank', {'Pass rate','Median error'}, 'Pass rate', 'How the best position is chosen.', 'none');
addButtons(gSw, {'Run sweep','Apply best','Export sweep CSV'}, ...
    {@(~,~)doSweep(), @(~,~)doApplyBestSweep(), @(~,~)doExportSweep()}, true);
lblSweep = addInfo(gSw, 'Re-runs the evaluation with the array moved to each candidate origin.', 2);

gDisp = mkPanelGrid(gFL, 'Display', 7, 3, 2, THEME.panelBg);
addDrop(gDisp, 'Colour by:', 'dispMetric', {'Score error (m)','Horizontal error (m)','3-D error (m)','Depth error (m)', ...
    'Bearing error (deg)','Min. SNR (dB)','# hydrophones detecting','Pass / fail','Sweep: pass rate (%)','Sweep: median error (m)'}, ...
    'Score error (m)', 'Quantity shown on the map.', 'disp');
addDrop(gDisp, 'View:', 'viewMode', {'Plan (map)','3-D'}, 'Plan (map)', 'Plan view or 3-D view with the seafloor surface.', 'disp');
addDrop(gDisp, 'Depth layer:', 'dispLayer', {'All depths'}, 'All depths', 'Which source-depth layer to show.', 'disp');
addNum(gDisp, 'Vertical exaggeration (3-D):', 'VE', 10, [0.1 1000], '%.1f', 'Vertical exaggeration of the 3-D view.', false, 'disp');
addNum(gDisp, 'Colour limit (0 = auto):', 'cmax', 0, [0 1e7], '%.1f', 'Upper colour limit; auto = 2 x threshold for errors.', false, 'disp');
addChk(gDisp, 'Show bathymetry', 'showBathy', true, 'Contours (plan) or seafloor surface (3-D).', 'disp');
addChk(gDisp, 'Label hydrophones', 'showLabels', true, 'Show H1..HN labels.', 'disp');

pFR = uipanel(gFoA, 'Title', 'FIELD OF ACCURACY', 'FontWeight', 'bold');
pFR.Layout.Row = 1; pFR.Layout.Column = 2; pFR.BackgroundColor = THEME.panelBg;
gFR = uigridlayout(pFR, [1 1]); gFR.Padding = [6 6 6 6];
axF = uiaxes(gFR);

%% -------------------- Initialise --------------------
applyPreset(true);
S.M = HydrophoneGeometry.pattern(CTRL.patName.Value, patternOpts());
refreshAll();
updateStatus('Ready. Design an array (tab 1), set the environment (tab 2) and calls (tab 3), then evaluate (tab 4).');

%% ====================================================================
%  UI helper constructors
%% ====================================================================
    function g = mkPanelGrid(parent, title, nRows, row, col, bg)
        p = uipanel(parent, 'Title', title, 'FontWeight', 'bold');
        p.Layout.Row = row; p.Layout.Column = col;
        p.BackgroundColor = bg;
        g = uigridlayout(p, [nRows 2]);
        g.RowHeight = repmat({22}, 1, nRows);
        g.ColumnWidth = {205, '1x'};
        g.Padding = [8 6 8 6];
        g.RowSpacing = 4;
        g.ColumnSpacing = 6;
        g.Scrollable = 'on';
        g.UserData = 0;
    end

    function r = nextRow(g, span)
        if nargin < 2, span = 1; end
        r0 = g.UserData + 1;
        g.UserData = g.UserData + span;
        if span > 1, r = [r0 r0+span-1]; else, r = r0; end
    end

    function h = addNum(g, label, key, def, lim, fmt, tip, isInt, grp)
        r = nextRow(g);
        lb = uilabel(g, 'Text', label, 'FontWeight', 'bold', 'Tooltip', tip);
        lb.Layout.Row = r; lb.Layout.Column = 1;
        h = uieditfield(g, 'numeric', 'Value', def, 'Limits', lim, 'ValueDisplayFormat', fmt, 'Tooltip', tip);
        if isInt, h.RoundFractionalValues = 'on'; end
        h.Layout.Row = r; h.Layout.Column = 2;
        h.Tag = key; h.ValueChangedFcn = PARAMCB;
        CTRL.(key) = h; DEF.(key) = def; GRP.(key) = grp;
    end

    function h = addDrop(g, label, key, items, def, tip, grp)
        r = nextRow(g);
        lb = uilabel(g, 'Text', label, 'FontWeight', 'bold', 'Tooltip', tip);
        lb.Layout.Row = r; lb.Layout.Column = 1;
        h = uidropdown(g, 'Items', items, 'Value', def, 'Tooltip', tip);
        h.Layout.Row = r; h.Layout.Column = 2;
        h.Tag = key; h.ValueChangedFcn = PARAMCB;
        CTRL.(key) = h; DEF.(key) = def; GRP.(key) = grp;
    end

    function h = addTxt(g, label, key, def, tip, grp)
        r = nextRow(g);
        lb = uilabel(g, 'Text', label, 'FontWeight', 'bold', 'Tooltip', tip);
        lb.Layout.Row = r; lb.Layout.Column = 1;
        h = uieditfield(g, 'text', 'Value', def, 'Tooltip', tip);
        h.Layout.Row = r; h.Layout.Column = 2;
        h.Tag = key; h.ValueChangedFcn = PARAMCB;
        CTRL.(key) = h; DEF.(key) = def; GRP.(key) = grp;
    end

    function h = addChk(g, label, key, def, tip, grp)
        r = nextRow(g);
        h = uicheckbox(g, 'Text', label, 'Value', def, 'Tooltip', tip);
        h.Layout.Row = r; h.Layout.Column = [1 2];
        h.Tag = key; h.ValueChangedFcn = PARAMCB;
        CTRL.(key) = h; DEF.(key) = def; GRP.(key) = grp;
    end

    function lb = addInfo(g, txt, span)
        r = nextRow(g, span);
        lb = uilabel(g, 'Text', txt, 'WordWrap', 'on', 'VerticalAlignment', 'top', 'FontColor', [0.25 0.25 0.3]);
        lb.Layout.Row = r; lb.Layout.Column = [1 2];
    end

    function gb = addButtons(g, labels, fcns, useRow)
        gb = uigridlayout(g, [1 numel(labels)]);
        gb.Padding = [0 0 0 0];
        gb.ColumnSpacing = 6;
        gb.ColumnWidth = repmat({'1x'}, 1, numel(labels));
        if useRow
            gb.Layout.Row = nextRow(g);
            gb.Layout.Column = [1 2];
        end
        for q = 1:numel(labels)
            b = uibutton(gb, 'Text', labels{q}, 'ButtonPushedFcn', fcns{q});
            b.BackgroundColor = THEME.secondaryBtn;
            b.FontColor = THEME.textColor;
        end
    end

    function b = mkBtn(g, txt, fcn, row, col, primary, tip)
        b = uibutton(g, 'Text', txt, 'ButtonPushedFcn', fcn, 'Tooltip', tip);
        b.Layout.Row = row; b.Layout.Column = col;
        if primary
            b.BackgroundColor = THEME.accentColor; b.FontColor = [1 1 1]; b.FontWeight = 'bold';
        else
            b.BackgroundColor = THEME.secondaryBtn; b.FontColor = THEME.textColor;
        end
    end

%% ====================================================================
%  Parameter handling
%% ====================================================================
    function onParamChanged(key)
        grp = GRP.(key);
        try
            switch grp
                case 'preset'
                    applyPreset(true);
                    clearResults();
                    refreshEnvPlots();
                case 'call'
                    updateCall(true);
                    clearResults();
                    refreshEnvPlots();
                case 'array'
                    clearResults();
                    refreshArrayTab(); refreshEnvPlots(); refreshFieldPlot();
                case 'env'
                    if any(strcmp(key, {'lat0','lon0'})), onOriginChanged(); end
                    clearResults();
                    refreshArrayTab(); refreshEnvPlots(); refreshGridPreview(); refreshFieldPlot();
                case 'field'
                    clearResults();
                    updateDepthLayers();
                    refreshGridPreview(); refreshFieldPlot(); refreshEnvPlots();
                case 'loc'
                    clearResults();
                    refreshEnvPlots(); refreshFieldPlot();
                case 'disp'
                    refreshFieldPlot();
                otherwise
                    % no immediate action
            end
            if any(strcmp(grp, {'call','array','env','field','loc','preset'}))
                updateStatus('Settings changed - evaluate the field again to update results.');
            end
        catch ME
            updateStatus(['Error: ' ME.message]);
        end
    end

    function P = readParams(doValidate)
        if nargin < 1, doValidate = true; end
        P = struct();
        % array
        P.ox = getv('ox'); P.oy = getv('oy');
        P.placeMode = placeCode(getv('placeMode'));
        P.yaw = getv('yaw'); P.pitch = getv('pitch'); P.roll = getv('roll');
        P.refDepth = getv('refDepth');
        P.Lmax = getv('Lmax'); P.minPair = getv('minPair');
        % environment
        P.flatDepth = getv('flatDepth');
        P.sspModel = getv('sspModel'); P.c0 = getv('c0');
        switch getv('locSSPMode')
            case 'Constant c (below)', P.locSSPMode = 'constant';
            case 'Ray-traced (profile at array)', P.locSSPMode = 'ray';
            otherwise, P.locSSPMode = 'env';
        end
        if startsWith(getv('propModel'), 'Ray'), P.propModel = 'ray'; else, P.propModel = 'straight'; end
        P.BL = getv('BL'); P.nRays = round(getv('nRays'));
        switch getv('rayArrivals')
            case 'Strongest arrival only', P.rayArrivals = 'strongest';
            case 'First arrival only', P.rayArrivals = 'first';
            otherwise, P.rayArrivals = 'all';
        end
        P.cLoc = getv('cLoc');
        if strcmp(getv('snrMode'), 'Fixed SNR'), P.snrMode = 'fixed'; else, P.snrMode = 'sonar'; end
        P.snrFixed = getv('snrFixed'); P.NL = getv('NL');
        sp = getv('spreading');
        if startsWith(sp, 'Spherical'), P.spreading = 'spherical';
        elseif startsWith(sp, 'Practical'), P.spreading = 'practical';
        elseif startsWith(sp, 'Cylindrical'), P.spreading = 'cylindrical';
        else, P.spreading = 'hybrid';
        end
        P.transR = getv('transR'); P.absorption = logical(getv('absorption'));
        P.DT = getv('DT'); P.minDetect = round(getv('minDetect'));
        P.timingSigma = getv('timingMs')/1000; P.posSigma = getv('posSigma');
        % call
        P.preset = getv('callPreset');
        P.fs = getv('fs'); P.dur = getv('dur'); P.f0 = getv('f0'); P.f1 = getv('f1');
        P.sweep = getv('sweep'); P.taper = getv('taper')/100; P.SL = getv('SL'); P.tplFs = getv('tplFs');
        % field
        P.fcx = getv('fcx'); P.fcy = getv('fcy'); P.R = getv('R'); P.Rin = getv('Rin'); P.step = getv('step');
        if strcmp(getv('layout'), 'Polar rings'), P.layout = 'polar'; else, P.layout = 'cartesian'; end
        P.azFrom = getv('azFrom'); P.azTo = getv('azTo');
        P.depths = parseDepths(getv('depths'));
        P.minDistAny = getv('minDistAny');
        % localisation
        if startsWith(getv('solveMode'), '3-D'), P.solveMode = '3d'; else, P.solveMode = '2d'; end
        P.assumedDepth = getv('assumedDepth');
        P.useTrueDepth = logical(getv('useTrueDepth'));
        if startsWith(getv('tdoaMethod'), 'Envelope'), P.tdoaMethod = 'envelope'; else, P.tdoaMethod = 'waveform'; end
        if startsWith(getv('scoreMetric'), '3-D'), P.scoreMetric = '3d'; else, P.scoreMetric = 'horizontal'; end
        if startsWith(getv('undetPolicy'), 'Exclude'), P.undetected = 'exclude'; else, P.undetected = 'fail'; end
        P.best = getv('best'); P.tol = getv('tol'); P.threshold_m = P.best + P.tol;
        P.lambda = getv('lambda'); P.failPen = getv('failPen'); P.seed = round(getv('seed'));
        P.useParallel = logical(getv('useParallel'));
        % optimiser
        P.maxIters = round(getv('maxIters')); P.probes = round(getv('probes'));
        P.initStep = getv('initStep'); P.minStep = getv('minStep'); P.stall = round(getv('stall'));
        P.shake = getv('shake'); P.epsJ = getv('epsJ'); P.targetPass = getv('targetPass');
        P.moveZ = startsWith(getv('moveMode'), '3-D');

        if doValidate && P.Rin >= P.R, error('Field inner radius must be smaller than the outer radius.'); end
        if strcmp(P.solveMode, '3d'), P.minDetect = max(P.minDetect, 4); else, P.minDetect = max(P.minDetect, 3); end
    end

    function val = getv(key)
        val = CTRL.(key).Value;
    end

    function a = parseAngles(str)
        parts = regexp(char(str), '[,;\s]+', 'split');
        a = str2double(parts);
        a = unique(wrapDeg(a(isfinite(a))));
    end

    function a = wrapDeg(a)
        a = mod(a + 180, 360) - 180;          % (-180, 180]
        a(a == -180) = 180;
    end

    function d = parseDepths(str)
        parts = regexp(char(str), '[,;\s]+', 'split');
        d = str2double(parts);
        d = abs(d(isfinite(d)));
        if isempty(d), d = 50; end
        d = unique(d);
    end

    function A = placement(P, ox, oy, yaw)
        if nargin < 2, ox = P.ox; oy = P.oy; end
        if nargin < 4, yaw = P.yaw; end
        A = struct('ox', ox, 'oy', oy, 'mode', P.placeMode, 'refDepth', P.refDepth, ...
            'yaw', yaw, 'pitch', P.pitch, 'roll', P.roll);
    end

    function c = placeCode(label)
        if startsWith(label, 'Seafloor frame'), c = 'frame';
        elseif startsWith(label, 'Fixed'), c = 'fixed';
        else, c = 'seafloor';
        end
    end

    function label = placeLabel(code)
        switch lower(code)
            case 'frame', label = 'Seafloor frame (rigid)';
            case 'fixed', label = 'Fixed depth (rigid)';
            otherwise, label = 'Seafloor-mounted (each element)';
        end
    end

    function env = currentEnv(P)
        env = S.env;
        env.flatDepth = P.flatDepth;
        switch P.sspModel
            case 'Constant'
                env.ssp = OceanEnvironment.constantSSP(P.c0);
            case '1-D profile (file)'
                if isempty(S.sspProfile), error('Truth model is "1-D profile" but no profile is loaded (Environment tab).'); end
                env.ssp = S.sspProfile;
            otherwise
                if isempty(S.ncom), error('Truth model is "3-D field" but no NCOM file is loaded (Environment tab).'); end
                env.ssp = OceanEnvironment.fieldSSP(S.ncom);
        end
        env.prop = struct('spreading', P.spreading, 'transR', P.transR, 'absorption', P.absorption);
    end

    function o = patternOpts()
        o = struct('N', CTRL.patN.Value, 'dx', CTRL.patDx.Value, 'dy', CTRL.patDy.Value, ...
            'nx', CTRL.patNx.Value, 'ny', CTRL.patNy.Value, 'rot', CTRL.patRot.Value, ...
            'jitter', CTRL.patJit.Value, 'z', CTRL.patZ.Value, 'Lmax', CTRL.Lmax.Value, 'minPair', CTRL.minPair.Value);
    end

    function ctx = buildContext(P)
        % Everything needed to evaluate a geometry: env, grid, cfg, scoring
        ctx = struct();
        ctx.P = P;
        ctx.env = currentEnv(P);
        if isempty(S.call), updateCall(false); end
        if isempty(S.call), error('No call signal available (Calls tab).'); end
        F = struct('cx', P.fcx, 'cy', P.fcy, 'R', P.R, 'Rin', P.Rin, 'step', P.step, ...
            'layout', P.layout, 'azFrom', P.azFrom, 'azTo', P.azTo, 'depths', P.depths);
        [ctx.grid, ctx.gridInfo] = wahi_build_field_grid(F, ctx.env);
        if isempty(ctx.grid), error('No synthetic sources: field is empty, on land, or deeper than the seabed.'); end
        cfg = struct();
        cfg.fs = P.fs; cfg.call = S.call; cfg.band = S.band;
        cfg.snrMode = P.snrMode; cfg.snrFixed = P.snrFixed; cfg.SL = P.SL; cfg.NL = P.NL;
        cfg.DT = P.DT; cfg.minDetect = P.minDetect;
        cfg.timingSigma = P.timingSigma; cfg.posSigma = P.posSigma;
        cfg.solveMode = P.solveMode; cfg.assumedDepth = P.assumedDepth;
        cfg.tdoaMethod = P.tdoaMethod; cfg.seed = P.seed;
        cfg.useTrueDepth = P.useTrueDepth;
        cfg.locSSP = OceanEnvironment.constantSSP(P.cLoc);
        ctx.cfg = cfg;
        ctx.E = struct('threshold_m', P.threshold_m, 'failPenalty_m', P.failPen, ...
            'minDistToAnyHyd', P.minDistAny, 'scoreMetric', P.scoreMetric, ...
            'undetected', P.undetected, 'seed', P.seed, 'useParallel', P.useParallel, ...
            'progressFcn', PROGRESSCB);
    end

    function [st, dg, Hn] = evaluateGeometry(Mrel, ox, oy, ctx, yaw)
        P = ctx.P;
        if nargin < 5, yaw = P.yaw; end
        Hn = HydrophoneGeometry.toAbsolute(Mrel, placement(P, ox, oy, yaw), ctx.env);
        wd = OceanEnvironment.waterDepth(ctx.env, Hn(:,1), Hn(:,2));
        if any(wd <= 0), error('A hydrophone is on land at this position.'); end
        if any(HydrophoneGeometry.belowSeabed(Hn, ctx.env)), error('A hydrophone is below the seabed (check pitch/roll / reference depth).'); end
        Xp = Hn(:,1:2) - mean(Hn(:,1:2), 1);
        if rank(Xp, 1e-6*max(1, max(abs(Xp(:))))) < 2, error('Hydrophones are collinear in plan view (no horizontal fix).'); end
        cfg = wahi_geometry_cfg(ctx.cfg, P, ctx.env, Hn, ctx.grid);
        loc = WhaleCallLocaliser(cfg, ctx.env, Hn);
        [st, dg] = wahi_evaluate_field(loc, ctx.grid, ctx.E);
    end

    function abort = progressCb(i, n)
        updateStatus(sprintf('%s source %d / %d', S.progressPrefix, i, n));
        drawnow limitrate;
        abort = S.stopFlag;
    end

    function J = scoreJ(st, P)
        J = st.passRate - P.lambda * min(st.meanErr / P.failPen, 1);
    end

%% ====================================================================
%  Calls
%% ====================================================================
    function applyPreset(doUpdate)
        p = WhaleCallLibrary.getPreset(CTRL.callPreset.Value);
        lblCallNote.Text = p.note;
        if ~startsWith(p.kind, 'template')
            CTRL.f0.Value = p.f0; CTRL.f1.Value = p.f1; CTRL.dur.Value = p.dur;
            CTRL.fs.Value = p.fs; CTRL.sweep.Value = p.sweep; CTRL.SL.Value = p.SL;
        elseif isempty(S.callTemplate.values)
            updateStatus('Template preset selected - click "Import template..." to load a file.');
        end
        if nargin < 1 || doUpdate, updateCall(true); else, updateCall(false); end
    end

    function setPresetNote()
        pp = WhaleCallLibrary.getPreset(CTRL.callPreset.Value);
        lblCallNote.Text = pp.note;
    end

    function spec = callSpec()
        p = WhaleCallLibrary.getPreset(CTRL.callPreset.Value);
        spec = struct('kind', p.kind, 'f0', CTRL.f0.Value, 'f1', CTRL.f1.Value, 'dur', CTRL.dur.Value, ...
            'sweep', CTRL.sweep.Value, 'taper', CTRL.taper.Value/100, 'template', [], 'templateFs', NaN, 'templateAmp', []);
        if startsWith(p.kind, 'template')
            spec.template = S.callTemplate.values;
            spec.templateAmp = S.callTemplate.amp;
            spec.templateFs = S.callTemplate.fs;
            if ~isfinite(spec.templateFs) && CTRL.tplFs.Value > 0, spec.templateFs = CTRL.tplFs.Value; end
        end
    end

    function updateCall(doPlot)
        if nargin < 1, doPlot = true; end
        spec = callSpec();
        fs = CTRL.fs.Value;
        if startsWith(spec.kind, 'template') && isempty(spec.template)
            S.call = []; S.band = [NaN NaN];
            lblBand.Text = 'Band: --- (no template loaded)';
            if doPlot, plotCall(); end
            return;
        end
        if ~startsWith(spec.kind, 'template') && max(spec.f0, spec.f1) >= fs/2
            updateStatus(sprintf('Warning: call frequency %.0f Hz is above Nyquist (fs/2 = %.0f Hz) - increase the sample rate.', max(spec.f0,spec.f1), fs/2));
        end
        S.call = WhaleCallLibrary.synthesize(spec, fs);
        [a, b] = WhaleCallLibrary.estimateBand(S.call, fs);
        S.band = [a b];
        lblBand.Text = sprintf('Band used for noise/filtering: %.4g - %.4g Hz  |  %d samples (%.4g s)', a, b, numel(S.call), numel(S.call)/fs);
        if doPlot, plotCall(); end
    end

    function plotCall()
        cla(axWave); cla(axSpec);
        if isempty(S.call)
            title(axWave, 'No call'); return;
        end
        fs = CTRL.fs.Value;
        t = (0:numel(S.call)-1)/fs;
        plot(axWave, t, S.call, 'Color', THEME.accentColor);
        xlabel(axWave, 'Time (s)'); ylabel(axWave, 'Amplitude');
        title(axWave, sprintf('%s', CTRL.callPreset.Value), 'Interpreter', 'none');
        grid(axWave, 'on');
        axis(axWave, 'tight');
        [Sp, f, tt] = WhaleCallLibrary.spectrogramSimple(S.call, fs);
        imagesc(axSpec, tt, f, Sp);
        axSpec.YDir = 'normal';
        colormap(axSpec, parula);
        mx = max(Sp(:));
        caxis(axSpec, [mx-60 mx]);
        fTop = min(fs/2, max(S.band(2)*1.5, S.band(2) + 10));
        ylim(axSpec, [0 fTop]);
        xlim(axSpec, [0 max(numel(S.call)/fs, tt(end))]);
        xlabel(axSpec, 'Time (s)'); ylabel(axSpec, 'Frequency (Hz)');
        title(axSpec, 'Spectrogram (dB, 60 dB range)');
    end

    function doImportTemplate()
        [file, pth] = uigetfile({'*.csv;*.txt;*.mat;*.wav','Template files (*.csv, *.txt, *.mat, *.wav)'; '*.*','All files'}, ...
            'Select call template');
        figure(fig);
        if isequal(file, 0), updateStatus('Template import cancelled'); return; end
        try
            tpl = WhaleCallLibrary.loadTemplate(fullfile(pth, file));
            if strcmp(tpl.guess, 'contour'), defOpt = 1; else, defOpt = 2; end
            if isfinite(tpl.fs) && endsWith(lower(file), '.wav')
                choice = 'Waveform samples';
            else
                msg = sprintf(['Loaded %d values (range %.4g to %.4g).\n\n' ...
                    'How should they be interpreted?\n' ...
                    '- Frequency contour: values are Hz, stretched over the Call duration.\n' ...
                    '- Waveform samples: pressure samples at the Template sample rate (%s).'], ...
                    numel(tpl.values), min(tpl.values), max(tpl.values), fsText());
                choice = uiconfirm(fig, msg, 'Template type', 'Options', {'Frequency contour','Waveform samples','Cancel'}, ...
                    'DefaultOption', defOpt, 'CancelOption', 3);
            end
            if strcmp(choice, 'Cancel'), updateStatus('Template import cancelled'); return; end
            S.callTemplate = tpl;
            if strcmp(choice, 'Frequency contour')
                CTRL.callPreset.Value = 'Imported template - freq. contour';
                if max(tpl.values) >= CTRL.fs.Value/2
                    CTRL.fs.Value = ceil(2.5*max(tpl.values));
                    updateStatus(sprintf('Sample rate raised to %.0f Hz to cover the contour.', CTRL.fs.Value));
                end
                lblTpl.Text = sprintf('Contour: %s (%d points, %.4g-%.4g Hz)', file, numel(tpl.values), min(tpl.values), max(tpl.values));
            else
                CTRL.callPreset.Value = 'Imported template - waveform';
                S.callTemplate.amp = [];
                if isfinite(tpl.fs)
                    lblTpl.Text = sprintf('Waveform: %s (%d samples at %.0f Hz)', file, numel(tpl.values), tpl.fs);
                else
                    lblTpl.Text = sprintf('Waveform: %s (%d samples at %s)', file, numel(tpl.values), fsText());
                end
            end
            setPresetNote();
            updateCall(true);
            clearResults();
            refreshEnvPlots();
            updateStatus(sprintf('Template loaded: %s', file));
        catch ME
            uialert(fig, ME.message, 'Template error');
        end
        function s = fsText()
            if CTRL.tplFs.Value > 0, s = sprintf('%.0f Hz', CTRL.tplFs.Value);
            else, s = sprintf('fs = %.0f Hz', CTRL.fs.Value);
            end
        end
    end

%% ====================================================================
%  Array design tab
%% ====================================================================
    function doGeneratePattern()
        try
            Mn = HydrophoneGeometry.pattern(CTRL.patName.Value, patternOpts());
            S.M = Mn;
            % relax the optimiser constraints so the new design is valid
            note = '';
            Dp = HydrophoneGeometry.pairDist(Mn); Dp(1:size(Dp,1)+1:end) = inf;
            rmaxP = max(sqrt(sum((Mn(:,1:2) - mean(Mn(:,1:2),1)).^2, 2)));
            if min(Dp(:)) < CTRL.minPair.Value
                CTRL.minPair.Value = max(0, floor(0.9*min(Dp(:))*1000)/1000);
                note = sprintf(' Min spacing set to %.3g m.', CTRL.minPair.Value);
            end
            if rmaxP > CTRL.Lmax.Value
                CTRL.Lmax.Value = min(1e5, max(1, ceil(1.25*rmaxP*1000)/1000));
                note = sprintf('%s Lmax set to %.3g m.', note, CTRL.Lmax.Value);
            end
            clearResults();
            refreshArrayTab(); refreshEnvPlots(); refreshFieldPlot();
            [ok, msg] = HydrophoneGeometry.isValid(S.M, CTRL.Lmax.Value, CTRL.minPair.Value);
            if ok, updateStatus(sprintf('Generated %d hydrophones (%s).%s', size(S.M,1), CTRL.patName.Value, note));
            else, updateStatus(['Generated pattern violates a constraint: ' msg]);
            end
        catch ME
            uialert(fig, ME.message, 'Pattern error');
        end
    end

    function onTableEdit(src, evt)
        D = src.Data;
        Dv = D(:,2:4);
        if any(~isfinite(Dv(:)))
            src.Data = tableData();
            updateStatus('Invalid value - edit reverted.');
            return;
        end
        S.M = D(:, 2:4);
        clearResults();
        refreshArrayTab(); refreshEnvPlots(); refreshFieldPlot();
        if ~isempty(evt) && isprop(evt, 'Indices') && ~isempty(evt.Indices)
            updateStatus(sprintf('Hydrophone %d updated.', evt.Indices(1)));
        end
    end

    function onTableSelect(~, evt)
        if isempty(evt.Indices), S.selRows = []; else, S.selRows = unique(evt.Indices(:,1))'; end
    end

    function D = tableData()
        n = size(S.M, 1);
        D = [(1:n)' S.M];
    end

    function doAddHyd()
        if isempty(S.M), S.M = zeros(0,3); end
        if size(S.M,1) >= 64, updateStatus('Maximum of 64 hydrophones.'); return; end
        if isempty(S.M)
            p = [0 0 CTRL.patZ.Value];
        else
            c = mean(S.M(:,1:2), 1);
            rmax = max(sqrt(sum((S.M(:,1:2) - c).^2, 2)));
            p = [c(1) + max(rmax, CTRL.minPair.Value)*0.5, c(2) + max(rmax, CTRL.minPair.Value)*0.5, S.M(end,3)];
        end
        S.M = [S.M; p];
        clearResults(); refreshArrayTab(); refreshFieldPlot();
        updateStatus(sprintf('Added H%d - edit its coordinates in the table.', size(S.M,1)));
    end

    function doRemoveHyd()
        if isempty(S.selRows), updateStatus('Select a cell in the row(s) to remove first.'); return; end
        if size(S.M,1) - numel(S.selRows) < 3, updateStatus('At least 3 hydrophones are required.'); return; end
        S.M(S.selRows, :) = [];
        S.selRows = [];
        clearResults(); refreshArrayTab(); refreshFieldPlot();
        updateStatus('Hydrophone(s) removed.');
    end

    function doRecentre()
        if isempty(S.M), return; end
        S.M = HydrophoneGeometry.recentreXY(S.M);
        clearResults(); refreshArrayTab(); refreshFieldPlot();
        updateStatus('Array re-centred on its horizontal centroid.');
    end

    function doRenumber()
        refreshArrayTab();
    end

    function doImportGeometry()
        [file, pth] = uigetfile({'*.csv;*.txt','Geometry CSV (*.csv, *.txt)'; '*.*','All files'}, 'Import hydrophone geometry');
        figure(fig);
        if isequal(file, 0), updateStatus('Import cancelled'); return; end
        try
            [Mn, meta] = HydrophoneGeometry.readCSV(fullfile(pth, file));
            if size(Mn,1) > 64, error('Maximum 64 hydrophones supported.'); end
            S.M = Mn;
            if isfinite(meta.ox) || ~isempty(meta.mode)
                choice = uiconfirm(fig, sprintf(['The file contains placement metadata (origin %.1f, %.1f m; mode "%s").\n' ...
                    'Apply it to the Placement settings?'], meta.ox, meta.oy, meta.mode), 'Apply placement?', ...
                    'Options', {'Apply','Keep current'}, 'DefaultOption', 1, 'CancelOption', 2);
                if strcmp(choice, 'Apply')
                    if isfinite(meta.ox), CTRL.ox.Value = meta.ox; end
                    if isfinite(meta.oy), CTRL.oy.Value = meta.oy; end
                    if ~isempty(meta.mode), CTRL.placeMode.Value = placeLabel(meta.mode); end
                    if isfinite(meta.yaw), CTRL.yaw.Value = wrapDeg(meta.yaw); end
                    if isfinite(meta.pitch), CTRL.pitch.Value = wrapDeg(meta.pitch); end
                    if isfinite(meta.roll), CTRL.roll.Value = wrapDeg(meta.roll); end
                    if isfinite(meta.refDepth), CTRL.refDepth.Value = meta.refDepth; end
                    if isfinite(meta.lat0) && isempty(strtrim(CTRL.lat0.Value))
                        CTRL.lat0.Value = sprintf('%.8f', meta.lat0); CTRL.lon0.Value = sprintf('%.8f', meta.lon0);
                        onOriginChanged();
                    end
                end
            end
            clearResults(); refreshAll();
            [ok, msg] = HydrophoneGeometry.isValid(S.M, CTRL.Lmax.Value, CTRL.minPair.Value);
            if ok, updateStatus(sprintf('Imported %d hydrophones from %s.', size(S.M,1), file));
            else, updateStatus(sprintf('Imported %d hydrophones - note: %s', size(S.M,1), msg));
            end
        catch ME
            uialert(fig, ME.message, 'Import error');
        end
    end

    function doExportGeometry()
        if isempty(S.M), uialert(fig, 'No geometry to export.', 'Nothing to export'); return; end
        [file, pth] = uiputfile({'*.csv','CSV (*.csv)'}, 'Export hydrophone geometry', ...
            sprintf('hydrophones_%dch_%s.csv', size(S.M,1), datestr(now, 'yyyymmdd_HHMMSS')));
        figure(fig);
        if isequal(file, 0), updateStatus('Export cancelled'); return; end
        try
            P = readParams();
            HydrophoneGeometry.writeCSV(fullfile(pth, file), S.M, placement(P), currentEnvSafe(P));
            updateStatus(sprintf('Saved %d hydrophones to %s (relative + absolute coordinates).', size(S.M,1), file));
        catch ME
            uialert(fig, ME.message, 'Export error');
        end
    end

    function env = currentEnvSafe(P)
        try
            env = currentEnv(P);
        catch
            env = S.env; env.flatDepth = P.flatDepth;
        end
    end

    function goEvaluate()
        tg.SelectedTab = tabFoA;
        doRunWAH();
    end

    function refreshArrayTab()
        if isempty(S.M), tblHyd.Data = zeros(0,4); return; end
        tblHyd.Data = tableData();
        P = readParamsSafe();
        switch P.placeMode
            case 'seafloor', tblHyd.ColumnName = {'#','x (m)','y (m)','z HAB (m)'};
            case 'frame',    tblHyd.ColumnName = {'#','x (m)','y (m)','z above base (m)'};
            otherwise,       tblHyd.ColumnName = {'#','x (m)','y (m)','z up (m)'};
        end
        env = currentEnvSafe(P);
        A = placement(P);
        Hn = HydrophoneGeometry.toAbsolute(S.M, A, env);
        W = [Hn(:,1) - P.ox, Hn(:,2) - P.oy];          % deployed plan offsets
        n = size(S.M,1);
        cmap = lines(n);

        % ---- plan view (after orientation) ----
        ax = axArrPlan;
        cla(ax); resetAxesModes(ax); hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
        c = mean(W, 1);
        dxy = abs(W - c);
        ext = max([max(dxy(:)), CTRL.minPair.Value, 0.5]);
        if cbArrLmax.Value, ext = max(ext, CTRL.Lmax.Value); end
        lim = 1.2*ext;
        if cbArrBathy.Value && ~isempty(env.bathy)
            [Xc, Yc, Dc] = bathyWindow(env.bathy, P.ox + c(1) + [-lim lim], P.oy + c(2) + [-lim lim]);
            if ~isempty(Dc) && numel(unique(round(Dc(:)))) > 1
                [cc, hc] = contour(ax, Xc - P.ox, Yc - P.oy, Dc, 8, 'LineColor', [0.7 0.7 0.75]);
                try, clabel(cc, hc, 'Color', [0.5 0.5 0.55], 'FontSize', 8); catch, end
            end
        end
        if cbArrLmax.Value
            a = linspace(0, 2*pi, 181);
            plot(ax, c(1) + CTRL.Lmax.Value*cos(a), c(2) + CTRL.Lmax.Value*sin(a), '--', 'Color', [0.6 0.6 0.6]);
        end
        Dm = HydrophoneGeometry.pairDist(Hn);
        if cbPairLines.Value && n <= 16
            for a1 = 1:n
                for a2 = a1+1:n
                    plot(ax, W([a1 a2],1), W([a1 a2],2), '-', 'Color', [0.75 0.82 0.9]);
                    if n <= 8
                        text(ax, mean(W([a1 a2],1)), mean(W([a1 a2],2)), fmtLen(Dm(a1,a2)), ...
                            'FontSize', 8, 'Color', [0.35 0.45 0.6], 'HorizontalAlignment', 'center', 'BackgroundColor', [1 1 1]);
                    end
                end
            end
        end
        scatter(ax, W(:,1), W(:,2), 90, cmap, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.2, 'Marker', 'v');
        for q = 1:n
            text(ax, W(q,1), W(q,2), sprintf('   H%d', q), 'FontWeight', 'bold', 'Color', cmap(q,:)*0.8);
        end
        plot(ax, 0, 0, 'p', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
        Rb = HydrophoneGeometry.orientMatrix(P.yaw, P.pitch, P.roll);
        if strcmp(P.placeMode, 'seafloor'), Rb = HydrophoneGeometry.orientMatrix(P.yaw, 0, 0); end
        ah = 0.35*ext;
        plot(ax, [0 Rb(1,1)]*ah, [0 Rb(2,1)]*ah, '-', 'Color', [0.85 0.2 0.15], 'LineWidth', 1.5);
        plot(ax, [0 Rb(1,2)]*ah, [0 Rb(2,2)]*ah, '-', 'Color', [0.1 0.6 0.2], 'LineWidth', 1.5);
        text(ax, Rb(1,1)*ah, Rb(2,1)*ah, ' body x', 'Color', [0.85 0.2 0.15], 'FontSize', 8);
        text(ax, Rb(1,2)*ah, Rb(2,2)*ah, ' body y', 'Color', [0.1 0.6 0.2], 'FontSize', 8);
        xlim(ax, c(1) + [-lim lim]); ylim(ax, c(2) + [-lim lim]);
        xlabel(ax, 'East of origin (m)'); ylabel(ax, 'North of origin (m)');
        title(ax, sprintf('Plan view as deployed (yaw %.1f, pitch %.1f, roll %.1f deg) - origin (%.0f, %.0f) m', ...
            P.yaw, P.pitch, P.roll, P.ox, P.oy));
        hold(ax, 'off');

        % ---- 3-D view (world offsets, elevation relative to seafloor at origin) ----
        ax = axArrSide;
        cla(ax); resetAxesModes(ax); hold(ax, 'on'); grid(ax, 'on');
        wd0 = OceanEnvironment.waterDepth(env, P.ox, P.oy);
        zRel = Hn(:,3) + wd0;                          % height above seafloor at the origin
        gx = c(1) + linspace(-lim, lim, 25); gy = c(2) + linspace(-lim, lim, 25);
        [GX, GY] = meshgrid(gx, gy);
        GZ = -OceanEnvironment.waterDepth(env, P.ox + GX, P.oy + GY) + wd0;
        surf(ax, GX, GY, GZ, 'FaceColor', [0.82 0.76 0.64], 'EdgeColor', 'none', 'FaceAlpha', 0.45);
        if n <= 12
            for a1 = 1:n
                for a2 = a1+1:n
                    plot3(ax, W([a1 a2],1), W([a1 a2],2), zRel([a1 a2]), '-', 'Color', [0.6 0.7 0.85]);
                end
            end
        end
        scatter3(ax, W(:,1), W(:,2), zRel, 70, cmap, 'filled', 'MarkerEdgeColor', 'k');
        for q = 1:n
            text(ax, W(q,1), W(q,2), zRel(q), sprintf('  H%d', q), 'FontSize', 8, 'Color', cmap(q,:)*0.8);
        end
        z0 = 0;                                        % body origin height above seafloor@origin
        if strcmp(P.placeMode, 'fixed'), z0 = wd0 - P.refDepth; end
        cols = {[0.85 0.2 0.15], [0.1 0.6 0.2], [0.15 0.35 0.85]};
        for q = 1:3
            plot3(ax, [0 Rb(1,q)]*ah, [0 Rb(2,q)]*ah, z0 + [0 Rb(3,q)]*ah, '-', 'Color', cols{q}, 'LineWidth', 2);
        end
        vSpan = max(zRel) - min([zRel; GZ(:)]);
        if vSpan > 0.02*(2*lim)
            axis(ax, 'equal'); ttl = '3-D view (true scale)';
        else
            ttl = '3-D view (vertical exaggerated)';
        end
        view(ax, 35, 25);
        xlabel(ax, 'East (m)'); ylabel(ax, 'North (m)'); zlabel(ax, 'Height above seafloor at origin (m)');
        title(ax, [ttl ' - red/green/blue = body x/y/z']);
        hold(ax, 'off');

        % ---- tables and stats ----
        tblWorld.Data = [(1:n)' round(W(:,1),3) round(W(:,2),3) round(zRel,3) round(-Hn(:,3),3)];
        tblWorld.ColumnName = {'#','East (m)','North (m)','Above seafloor@origin (m)','Depth (m)'};
        tblPair.Data = round(Dm, 3);
        tblPair.ColumnName = compose('H%d', 1:n);
        tblPair.RowName = compose('H%d', 1:n);
        Dn = Dm; Dn(1:n+1:end) = inf;
        rr = sqrt(sum((W - c).^2, 2));
        [ok, msg] = HydrophoneGeometry.isValid(S.M, CTRL.Lmax.Value, CTRL.minPair.Value);
        bad = HydrophoneGeometry.belowSeabed(Hn, env);
        if any(bad)
            lblArrValid.Text = sprintf('! H%s below the seabed (check pitch/roll)', strjoin(compose('%d', find(bad)'), ', H'));
            lblArrValid.FontColor = [0.8 0.2 0.1];
        elseif ok
            lblArrValid.Text = ['Constraints OK ' msg]; lblArrValid.FontColor = [0.1 0.55 0.2];
        else
            lblArrValid.Text = ['! ' msg]; lblArrValid.FontColor = [0.8 0.2 0.1];
        end
        lblArrStats.Text = sprintf(['%d hydrophones  |  aperture (max 3-D pair) %s  |  min spacing %s\n' ...
            'max horizontal radius from centroid %s  |  depths %.2f - %.2f m'], n, fmtLen(max(Dm(:))), ...
            fmtLen(min(Dn(:))), fmtLen(max(rr)), min(-Hn(:,3)), max(-Hn(:,3)));
        if wd0 <= 0
            lblArrDepth.Text = sprintf('Origin is on LAND (elevation %.0f m) - move the array.', -wd0);
        else
            lblArrDepth.Text = sprintf('Water depth at origin: %.1f m', wd0);
        end
    end

    function s = fmtLen(v)
        if v < 10, s = sprintf('%.3f m', v); elseif v < 1000, s = sprintf('%.1f m', v); else, s = sprintf('%.0f m', v); end
    end

    function P = readParamsSafe()
        P = readParams(false);
    end

%% ====================================================================
%  Environment tab
%% ====================================================================
    function doLoadBathy()
        [file, pth] = uigetfile({'*.nc;*.grd;*.nc4;*.csv;*.txt;*.xyz;*.dat;*.mat', ...
            'Bathymetry (*.nc, *.grd, *.csv, *.txt, *.xyz, *.mat)'; '*.*','All files'}, 'Load bathymetry');
        figure(fig);
        if isequal(file, 0), updateStatus('Bathymetry load cancelled'); return; end
        try
            updateStatus('Loading bathymetry...'); drawnow;
            [lon0, lat0] = originLonLat();
            B = OceanEnvironment.loadBathymetry(fullfile(pth, file), lon0, lat0);
            S.env.bathy = B;
            if B.isGeo && ~(isfinite(lat0) && isfinite(lon0))
                CTRL.lat0.Value = sprintf('%.6f', B.lat0);
                CTRL.lon0.Value = sprintf('%.6f', B.lon0);
            end
            if B.isGeo, S.env.lat0 = B.lat0; S.env.lon0 = B.lon0; end
            updateBathyInfo();
            clearResults(); refreshAll();
            updateStatus(sprintf('Bathymetry loaded: %s', file));
        catch ME
            uialert(fig, ME.message, 'Bathymetry error');
            updateStatus('Bathymetry load failed');
        end
    end

    function updateBathyInfo()
        B = S.env.bathy;
        if isempty(B)
            lblBathy.Text = 'No bathymetry loaded - flat bottom.';
            return;
        end
        [~, nm, ex] = fileparts(B.file);
        if B.isGeo, geo = sprintf('geographic, origin %.4f N %.4f E', B.lat0, B.lon0); else, geo = 'local metres'; end
        if B.signFlipped, sgn = 'elevation -> depth'; else, sgn = 'depth'; end
        lblBathy.Text = sprintf('%s%s (%s; %s)\nx %.1f..%.1f km, y %.1f..%.1f km, depth %.0f..%.0f m', nm, ex, geo, sgn, ...
            min(B.x)/1e3, max(B.x)/1e3, min(B.y)/1e3, max(B.y)/1e3, B.minDepth, B.maxDepth);
    end

    function doFlipBathy()
        if isempty(S.env.bathy), updateStatus('No bathymetry loaded.'); return; end
        S.env.bathy = OceanEnvironment.flipBathySign(S.env.bathy);
        updateBathyInfo(); clearResults(); refreshAll();
        updateStatus('Bathymetry depth sign flipped.');
    end

    function doClearBathy()
        S.env.bathy = [];
        updateBathyInfo(); clearResults(); refreshAll();
        updateStatus('Bathymetry cleared - using flat bottom.');
    end

    function [lon0, lat0] = originLonLat()
        lat0 = str2double(CTRL.lat0.Value); lon0 = str2double(CTRL.lon0.Value);
        if ~isfinite(lat0) || ~isfinite(lon0), lat0 = NaN; lon0 = NaN; end
    end

    function onOriginChanged()
        [lon0, lat0] = originLonLat();
        if ~isfinite(lat0)
            if ~isempty(S.env.bathy) && S.env.bathy.isGeo
                S.env.lat0 = S.env.bathy.lat0; S.env.lon0 = S.env.bathy.lon0;
                CTRL.lat0.Value = sprintf('%.6f', S.env.lat0); CTRL.lon0.Value = sprintf('%.6f', S.env.lon0);
            else
                S.env.lat0 = NaN; S.env.lon0 = NaN;
            end
            return;
        end
        S.env.lat0 = lat0; S.env.lon0 = lon0;
        try
            if ~isempty(S.env.bathy) && S.env.bathy.isGeo
                S.env.bathy = OceanEnvironment.loadBathymetry(S.env.bathy.file, lon0, lat0);
                updateBathyInfo();
            end
            if ~isempty(S.ncom)
                S.ncom = OceanEnvironment.loadNCOM(S.ncom.file, lon0, lat0, round(CTRL.ncomT.Value));
            end
            updateStatus('Geographic origin changed - grids re-projected.');
        catch ME
            updateStatus(['Re-projection failed: ' ME.message]);
        end
    end

    function doLoadProfile()
        [file, pth] = uigetfile({'*.csv;*.txt;*.dat;*.mat','Sound speed profile (*.csv, *.txt, *.mat)'; '*.*','All files'}, ...
            'Load 1-D sound speed profile');
        figure(fig);
        if isequal(file, 0), updateStatus('Profile load cancelled'); return; end
        try
            S.sspProfile = OceanEnvironment.loadProfile(fullfile(pth, file));
            CTRL.sspModel.Value = '1-D profile (file)';
            lblSSP.Text = sprintf('Profile %s: %d levels to %.0f m, c %.1f-%.1f m/s (%s)', file, numel(S.sspProfile.z), ...
                max(S.sspProfile.z), min(S.sspProfile.c), max(S.sspProfile.c), S.sspProfile.note);
            clearResults(); refreshEnvPlots();
            updateStatus(sprintf('Sound speed profile loaded: %s', file));
        catch ME
            uialert(fig, ME.message, 'Profile error');
        end
    end

    function doLoadNCOM()
        [file, pth] = uigetfile({'*.nc;*.nc4;*.cdf','NetCDF (*.nc, *.nc4)'; '*.*','All files'}, 'Load NCOM / HYCOM NetCDF');
        figure(fig);
        if isequal(file, 0), updateStatus('NCOM load cancelled'); return; end
        try
            updateStatus('Loading NCOM field (this can take a moment)...'); drawnow;
            [lon0, lat0] = originLonLat();
            S.ncom = OceanEnvironment.loadNCOM(fullfile(pth, file), lon0, lat0, round(CTRL.ncomT.Value));
            if ~isfinite(lat0)
                CTRL.lat0.Value = sprintf('%.6f', S.ncom.lat0); CTRL.lon0.Value = sprintf('%.6f', S.ncom.lon0);
                S.env.lat0 = S.ncom.lat0; S.env.lon0 = S.ncom.lon0;
            end
            CTRL.sspModel.Value = '3-D field (NCOM NetCDF)';
            lblSSP.Text = sprintf('NCOM %s: %d x %d profiles, %d levels to %.0f m; %s', file, numel(S.ncom.lon), ...
                numel(S.ncom.lat), numel(S.ncom.z), max(S.ncom.z), S.ncom.note);
            clearResults(); refreshEnvPlots();
            updateStatus(sprintf('NCOM field loaded: %s', file));
        catch ME
            uialert(fig, ME.message, 'NCOM error');
            updateStatus('NCOM load failed');
        end
    end

    function refreshEnvPlots()
        P = readParamsSafe();
        env = currentEnvSafe(P);
        % --- bathymetry map ---
        ax = axBathy;
        cla(ax); colorbar(ax, 'off');
        resetAxesModes(ax);
        hold(ax, 'on');
        if ~isempty(env.bathy)
            B = env.bathy;
            [xd, yd, Dd] = decimateForPlot(B, 400);
            imagesc(ax, xd/1e3, yd/1e3, Dd);
            ax.YDir = 'normal';
            colormap(ax, flipud(parula));
            cb = colorbar(ax); cb.Label.String = 'Water depth (m)';
            try
                contour(ax, xd/1e3, yd/1e3, Dd, [0 0], 'LineColor', 'k', 'LineWidth', 1.5);
            catch
            end
            ttl = 'Bathymetry';
        else
            ttl = sprintf('Flat bottom (%.0f m)', env.flatDepth);
        end
        drawFieldOutline(ax, P, 'k');
        if ~isempty(S.M)
            Hn = HydrophoneGeometry.toAbsolute(S.M, placement(P), env);
            plot(ax, Hn(:,1)/1e3, Hn(:,2)/1e3, 'v', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
        end
        axis(ax, 'equal'); grid(ax, 'on');
        xlabel(ax, 'East (km)'); ylabel(ax, 'North (km)'); title(ax, ttl);
        if isempty(env.bathy)
            Rv = max(CTRL.R.Value, 1000)/1e3*1.1;
            xlim(ax, CTRL.fcx.Value/1e3 + [-Rv Rv]); ylim(ax, CTRL.fcy.Value/1e3 + [-Rv Rv]);
        end
        hold(ax, 'off');

        % --- sound speed profiles ---
        ax = axSSP;
        cla(ax); hold(ax, 'on'); grid(ax, 'on');
        dmax = OceanEnvironment.maxWaterDepth(env);
        hL = gobjects(0); leg = {};
        if strcmp(env.ssp.type, 'field')
            fld = env.ssp.field;
            nx = numel(fld.x); ny = numel(fld.y);
            pick = unique(round(linspace(1, nx*ny, min(25, nx*ny))));
            for q = pick
                [ix, iy] = ind2sub([nx ny], q);
                plot(ax, squeeze(fld.C(ix, iy, :)), fld.z, '-', 'Color', [0.75 0.75 0.8]);
            end
            pr = OceanEnvironment.profileAt(env, P.ox, P.oy);
            hL(end+1) = plot(ax, pr.c, pr.z, 'b-', 'LineWidth', 2);
            leg{end+1} = 'truth at array';
        else
            zz = linspace(0, dmax, 200);
            cc = interp1(env.ssp.Sz, env.ssp.cz, min(zz, env.ssp.Sz(end)));
            hL(end+1) = plot(ax, cc, zz, 'b-', 'LineWidth', 2);
            leg{end+1} = 'truth';
        end
        try
            if strcmp(P.locSSPMode, 'constant')
                hL(end+1) = plot(ax, [P.cLoc P.cLoc], [0 dmax], 'r--', 'LineWidth', 1.5);
            else
                pr = OceanEnvironment.profileAt(env, P.ox, P.oy);
                hL(end+1) = plot(ax, pr.c, pr.z, 'r--', 'LineWidth', 1.5);
            end
            leg{end+1} = 'localiser';
        catch
        end
        ax.YDir = 'reverse';
        ylim(ax, [0 max(dmax, 10)]);
        xlabel(ax, 'Sound speed (m/s)'); ylabel(ax, 'Depth (m)');
        title(ax, 'Sound speed');
        legend(ax, hL, leg, 'Location', 'southwest');
        hold(ax, 'off');

        % --- SNR vs range (and ray fan) ---
        ax = axSNR;
        cla(ax); legend(ax, 'off'); hold(ax, 'on'); grid(ax, 'on');
        rmax = max(1.2*(CTRL.R.Value + hypot(CTRL.fcx.Value - CTRL.ox.Value, CTRL.fcy.Value - CTRL.oy.Value)), 2000);
        rr = logspace(1, log10(rmax), 300);
        fc = mean(S.band);
        if ~isfinite(fc), fc = mean([CTRL.f0.Value CTRL.f1.Value]); end
        dList = parseDepths(CTRL.depths.Value);
        RP = rayPreview(P, env, rmax);
        hS = gobjects(0); labS = {};
        if strcmp(P.snrMode, 'fixed')
            snr = P.snrFixed*ones(size(rr));
        else
            snr = CTRL.SL.Value - OceanEnvironment.transmissionLoss(env, rr, fc) - CTRL.NL.Value;
        end
        isRay = strcmp(P.propModel, 'ray') && ~isempty(RP);
        if isRay
            hS(end+1) = plot(ax, rr/1e3, snr, '--', 'Color', [0.6 0.6 0.6]); labS{end+1} = 'straight-ray spreading law';
            try
                Sp = WhaleCallLocaliser.cachedPaths(RP.ssp, RP.D, RP.zr, dList(1), RP.o);
                rl = linspace(10, rmax, 600);
                TLr = inf(size(rl));
                for q = 1:numel(rl)
                    Aq = OceanRayModel.eigenrays(Sp, rl(q), fc, P.absorption);
                    if ~isempty(Aq.t), TLr(q) = -10*log10(sum(10.^(-Aq.TL/10))); end
                end
                if strcmp(P.snrMode, 'fixed')
                    snrR = P.snrFixed*ones(size(rl)); snrR(~isfinite(TLr)) = NaN;
                else
                    snrR = CTRL.SL.Value - TLr - CTRL.NL.Value;
                end
                hS(end+1) = plot(ax, rl/1e3, snrR, 'b-', 'LineWidth', 1.5);
                labS{end+1} = sprintf('ray traced (H1 at %.0f m, source at %.0f m)', RP.zr, dList(1));
                rr = rl; snr = snrR;
            catch ME
                title(ax, ['Ray SNR failed: ' ME.message]);
            end
        else
            hS(end+1) = plot(ax, rr/1e3, snr, 'b-', 'LineWidth', 2); labS{end+1} = 'SNR';
        end
        hS(end+1) = plot(ax, [10 rmax]/1e3, CTRL.DT.Value*[1 1], 'r--'); labS{end+1} = 'detection threshold';
        if isRay, ax.XScale = 'linear'; else, ax.XScale = 'log'; end
        xlabel(ax, 'Horizontal range (km)'); ylabel(ax, 'Received SNR (dB)');
        det = snr >= CTRL.DT.Value;
        iDet = find(det, 1, 'last');
        if isempty(iDet), dr = 'never detected';
        elseif iDet == numel(rr), dr = sprintf('detected beyond %.0f km', rr(end)/1e3);
        else, dr = sprintf('last detection ~%.1f km', rr(iDet)/1e3);
        end
        if isRay && any(~det(1:max(iDet,1))), dr = [dr ', with gaps (shadow zones)']; end
        title(ax, sprintf('SNR vs range at %.4g Hz (SL %.0f, NL %.0f dB): %s', fc, CTRL.SL.Value, CTRL.NL.Value, dr));
        legend(ax, hS, labS, 'Location', 'northeast');
        hold(ax, 'off');

        ax = axRays;
        cla(ax); hold(ax, 'on'); grid(ax, 'on');
        if isempty(RP)
            title(ax, 'Ray paths: no hydrophones');
        else
            try
                Fr = OceanRayModel.fan(RP.ssp, RP.D, RP.zr, RP.o);
                ang = [0.5 1 2 3 5 8 12 17 25 35];
                rmaxP = min(rmax, 100e3);
                Pp = OceanRayModel.rayPaths(Fr, ang, rmaxP);
                cm = parula(numel(Pp) + 2);
                for q = 1:numel(Pp)
                    plot(ax, Pp{q}(:,1)/1e3, Pp{q}(:,2), '-', 'Color', cm(q,:), 'LineWidth', 0.8);
                end
                plot(ax, [0 rmaxP]/1e3, RP.D*[1 1], '-', 'Color', [0.45 0.35 0.2], 'LineWidth', 2);
                for q = 1:numel(dList)
                    plot(ax, [0 rmaxP]/1e3, dList(q)*[1 1], ':', 'Color', [0.8 0.2 0.1]);
                end
                ax.YDir = 'reverse';
                xlim(ax, [0 rmaxP/1e3]); ylim(ax, [0 RP.D*1.02]);
                xlabel(ax, 'Range from H1 (km)'); ylabel(ax, 'Depth (m)');
                title(ax, sprintf('Ray paths from H1 (%.0f m) at %s deg, up & down - profile at array, flat bottom %.0f m; dotted = source depths', ...
                    RP.zr, strjoin(compose('%g', ang), '/'), RP.D));
            catch ME
                title(ax, ['Ray plot failed: ' ME.message]);
            end
        end
        hold(ax, 'off');
    end

    function RP = rayPreview(P, env, rmax)
        % Ray-model settings for the plots (same rules as the evaluation)
        RP = [];
        if isempty(S.M), return; end
        try
            Hn = HydrophoneGeometry.toAbsolute(S.M, placement(P), env);
            dList = parseDepths(CTRL.depths.Value);
            G = [P.fcx + [0; rmax*0.8] P.fcy*[1; 1] -max(dList)*[1; 1]];
            cfgP = wahi_geometry_cfg(struct(), P, env, Hn, G);
            RP = cfgP.locRay;
            RP.o.rmax = 10000*ceil(rmax/10000);
            RP.zr = -Hn(1,3);
        catch
            RP = [];
        end
    end

%% ====================================================================
%  Field tab: grid preview, plot
%% ====================================================================
    function updateDepthLayers()
        d = parseDepths(CTRL.depths.Value);
        items = [{'All depths'}, compose('%g m', d)];
        cur = CTRL.dispLayer.Value;
        CTRL.dispLayer.Items = items;
        if any(strcmp(items, cur)), CTRL.dispLayer.Value = cur; else, CTRL.dispLayer.Value = 'All depths'; end
    end

    function refreshGridPreview()
        try
            P = readParams();
            env = currentEnv(P);
            F = struct('cx', P.fcx, 'cy', P.fcy, 'R', P.R, 'Rin', P.Rin, 'step', P.step, ...
                'layout', P.layout, 'azFrom', P.azFrom, 'azTo', P.azTo, 'depths', P.depths);
            [S.gridPts, S.gridInfo] = wahi_build_field_grid(F, env);
            txt = sprintf('%d sources (%s per depth)', size(S.gridPts,1), strjoin(compose('%d', S.gridInfo.nPerDepth), ' / '));
            if S.gridInfo.nDroppedSeabed > 0
                txt = sprintf('%s; %d dropped (land / below seabed)', txt, S.gridInfo.nDroppedSeabed);
            end
            if size(S.gridPts,1) > 3000
                txt = [txt '. Large grid - consider a coarser spacing for optimisation.'];
            end
            lblGrid.Text = txt;
        catch ME
            S.gridPts = []; S.gridInfo = [];
            lblGrid.Text = ['Grid error: ' ME.message];
        end
    end

    function [vals, lab, kind] = metricValues(st, key)
        kind = 'err';
        switch key
            case 'Score error (m)',        vals = st.score_m;  lab = 'Score error (m)';
            case 'Horizontal error (m)',   vals = st.errH_m;   lab = 'Horizontal error (m)';
            case '3-D error (m)',          vals = st.err3_m;   lab = '3-D error (m)';
            case 'Depth error (m)',        vals = st.errZ_m;   lab = 'Depth error (m)';
            case 'Bearing error (deg)',    vals = st.brg_deg;  lab = 'Bearing error (deg)'; kind = 'deg';
            case 'Min. SNR (dB)',          vals = st.minSNR;   lab = 'Min. received SNR (dB)'; kind = 'snr';
            case '# hydrophones detecting', vals = st.nDet;    lab = '# hydrophones detecting'; kind = 'count';
            otherwise,                     vals = st.score_m;  lab = 'Score error (m)';
        end
    end

    function refreshFieldPlot()
        ax = axF;
        cla(ax); colorbar(ax, 'off'); legend(ax, 'off');
        resetAxesModes(ax);
        hold(ax, 'on'); grid(ax, 'on');
        P = readParamsSafe();
        env = currentEnvSafe(P);
        is3 = strcmp(CTRL.viewMode.Value, '3-D');
        dispKey = CTRL.dispMetric.Value;

        % bathymetry
        if CTRL.showBathy.Value && ~isempty(env.bathy)
            [xd, yd, Dd] = decimateForPlot(env.bathy, 160);
            if is3
                surf(ax, xd/1e3, yd/1e3, -Dd, 'FaceColor', [0.82 0.76 0.64], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
            else
                if numel(unique(round(Dd(:)))) > 1
                    [cc, hc] = contour(ax, xd/1e3, yd/1e3, Dd, 10, 'LineColor', [0.72 0.72 0.76]);
                    try, clabel(cc, hc, 'Color', [0.55 0.55 0.6], 'FontSize', 7, 'LabelSpacing', 400); catch, end
                end
                try, contour(ax, xd/1e3, yd/1e3, Dd, [0 0], 'LineColor', [0.3 0.25 0.2], 'LineWidth', 1.5); catch, end
            end
        elseif is3 && CTRL.showBathy.Value
            Rv = (CTRL.R.Value + 1000)/1e3;
            [Xf, Yf] = meshgrid(CTRL.fcx.Value/1e3 + [-Rv Rv], CTRL.fcy.Value/1e3 + [-Rv Rv]);
            surf(ax, Xf, Yf, -env.flatDepth*ones(2), 'FaceColor', [0.82 0.76 0.64], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        end

        if ~is3, drawFieldOutline(ax, P, [0.35 0.35 0.35]); end

        % sources
        layerMask = [];
        if ~isempty(S.gridPts) && ~isempty(S.gridInfo)
            layerMask = true(size(S.gridPts,1), 1);
            if ~strcmp(CTRL.dispLayer.Value, 'All depths')
                dsel = str2double(strtok(CTRL.dispLayer.Value));
                k = find(abs(S.gridInfo.depths - dsel) < 1e-9, 1);
                if ~isempty(k), layerMask = S.gridInfo.depthIdx == k; end
            end
        end

        cbLabel = '';
        if startsWith(dispKey, 'Sweep')
            plotSweep(ax, dispKey, is3);
            cbLabel = 'sweep';
        elseif ~isempty(S.stats) && ~isempty(S.gridPts) && numel(S.stats.flag) == size(S.gridPts,1)
            st = S.stats;
            G = S.gridPts;
            used = st.usedMask & layerMask;
            failed = layerMask & (st.flag >= 2);
            if strcmp(dispKey, 'Pass / fail')
                pm = used & st.passMask; fm = used & ~st.passMask;
                scatterPts(ax, G(pm,:), 30, [0.15 0.7 0.3], is3);
                scatterPts(ax, G(fm,:), 30, [0.85 0.2 0.15], is3);
            else
                [vals, lab, kind] = metricValues(st, dispKey);
                show = used & isfinite(vals);
                if strcmp(kind, 'err') && ~strcmp(dispKey, 'Score error (m)'), show = show & st.flag == 1; end
                if any(show)
                    v = vals(show);
                    if CTRL.cmax.Value > 0
                        cl = [0 CTRL.cmax.Value];
                    else
                        switch kind
                            case 'err',   cl = [0 2*P.threshold_m];
                            case 'deg',   cl = [0 max(HydrophoneGeometry.prctileLocal(v, 95), 1e-3)];
                            case 'snr',   cl = [min(v) max(max(v), min(v)+1)];
                            otherwise,    cl = [0 max(max(v), 1)];
                        end
                    end
                    if strcmp(kind, 'snr') || strcmp(kind, 'count'), vv = v; else, vv = min(v, cl(2)); end
                    if is3
                        scatter3(ax, G(show,1)/1e3, G(show,2)/1e3, G(show,3), 30, vv, 'filled', 'MarkerFaceAlpha', 0.85);
                    else
                        scatter(ax, G(show,1)/1e3, G(show,2)/1e3, 34, vv, 'filled', 'MarkerFaceAlpha', 0.85);
                    end
                    colormap(ax, turbo);
                    caxis(ax, cl);
                    cbLabel = lab;
                end
            end
            if any(failed) && ~strcmp(dispKey, '# hydrophones detecting') && ~strcmp(dispKey, 'Min. SNR (dB)')
                scatterPts(ax, G(failed,:), 22, [0.55 0.55 0.55], is3, 'x');
            end
        elseif ~isempty(S.gridPts)
            G = S.gridPts(layerMask,:);
            scatterPts(ax, G, 10, [0.55 0.6 0.7], is3);
        end
        if ~isempty(cbLabel) && ~strcmp(cbLabel, 'sweep')
            cb = colorbar(ax); cb.Label.String = cbLabel;
        end

        % hydrophones
        if ~isempty(S.M)
            Hn = HydrophoneGeometry.toAbsolute(S.M, placement(P), env);
            if is3
                scatter3(ax, Hn(:,1)/1e3, Hn(:,2)/1e3, Hn(:,3), 70, 'k', 'filled', 'Marker', 'v');
            else
                scatter(ax, Hn(:,1)/1e3, Hn(:,2)/1e3, 60, 'k', 'filled', 'Marker', 'v');
            end
            if CTRL.showLabels.Value && size(Hn,1) <= 24
                for q = 1:size(Hn,1)
                    if is3
                        text(ax, Hn(q,1)/1e3, Hn(q,2)/1e3, Hn(q,3), sprintf(' H%d', q), 'FontSize', 8, 'FontWeight', 'bold');
                    else
                        text(ax, Hn(q,1)/1e3, Hn(q,2)/1e3, sprintf(' H%d', q), 'FontSize', 8, 'FontWeight', 'bold');
                    end
                end
            end
        end

        xlabel(ax, 'East (km)'); ylabel(ax, 'North (km)');
        if is3
            zlabel(ax, 'Elevation (m)');
            view(ax, -35, 28);
            daspect(ax, [1 1 1000/max(CTRL.VE.Value, 0.1)]);
        else
            view(ax, 2);
            axis(ax, 'equal');
        end
        title(ax, fieldTitle());
        hold(ax, 'off');
        drawnow limitrate;
    end

    function scatterPts(ax, G, sz, col, is3, mk)
        if isempty(G), return; end
        if nargin < 6, mk = 'o'; end
        if is3
            if strcmp(mk, 'x'), scatter3(ax, G(:,1)/1e3, G(:,2)/1e3, G(:,3), sz, col, 'Marker', 'x');
            else, scatter3(ax, G(:,1)/1e3, G(:,2)/1e3, G(:,3), sz, col, 'filled');
            end
        else
            if strcmp(mk, 'x'), scatter(ax, G(:,1)/1e3, G(:,2)/1e3, sz, col, 'Marker', 'x');
            else, scatter(ax, G(:,1)/1e3, G(:,2)/1e3, sz, col, 'filled');
            end
        end
    end

    function t = fieldTitle()
        if startsWith(CTRL.dispMetric.Value, 'Sweep')
            if isempty(S.sweep), t = 'No sweep yet - click "Run sweep"'; else, t = 'Array position sweep (each marker = array origin)'; end
        elseif isempty(S.stats)
            t = sprintf('Preview: %d synthetic sources (evaluate to colour by accuracy)', size(S.gridPts,1));
        else
            t = sprintf('%s  |  pass %.1f%% (threshold %.0f m)', CTRL.dispMetric.Value, 100*S.stats.passRate, S.stats.threshold_m);
        end
    end

    function plotSweep(ax, key, is3)
        if isempty(S.sweep), return; end
        sw = S.sweep;
        byMedian = strcmp(key, 'Sweep: median error (m)');
        % best candidate (over yaw) at each distinct position
        [U, ~, ic] = unique([sw.ox sw.oy], 'rows');
        v = nan(size(U,1),1); okU = false(size(U,1),1); yU = nan(size(U,1),1);
        for q = 1:size(U,1)
            m = find(ic == q & sw.valid);
            if isempty(m), continue; end
            if byMedian, [v(q), k] = min(sw.median(m)); else, [v(q), k] = max(100*sw.pass(m)); end
            yU(q) = sw.yaw(m(k)); okU(q) = true;
        end
        if byMedian, lab = 'Median error (m) with array here (best yaw)';
        else, lab = 'Pass rate (%) with array here (best yaw)';
        end
        if is3
            scatter3(ax, U(okU,1)/1e3, U(okU,2)/1e3, zeros(sum(okU),1), 120, v(okU), 'filled', 'Marker', 's');
        else
            scatter(ax, U(okU,1)/1e3, U(okU,2)/1e3, 140, v(okU), 'filled', 'Marker', 's', 'MarkerEdgeColor', [0.3 0.3 0.3]);
            if any(~okU), plot(ax, U(~okU,1)/1e3, U(~okU,2)/1e3, 'x', 'Color', [0.5 0.5 0.5]); end
            if numel(unique(sw.yaw)) > 1 && sum(okU) <= 60
                for q = find(okU)'
                    text(ax, U(q,1)/1e3, U(q,2)/1e3, sprintf('  %g deg', yU(q)), 'FontSize', 7);
                end
            end
            if isfinite(sw.best)
                plot(ax, sw.ox(sw.best)/1e3, sw.oy(sw.best)/1e3, 'p', 'MarkerSize', 18, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k');
            end
        end
        colormap(ax, turbo);
        if any(okU) && max(v(okU)) > min(v(okU)), caxis(ax, [min(v(okU)) max(v(okU))]); end
        cb = colorbar(ax); cb.Label.String = lab;
    end

    function drawFieldOutline(ax, P, col)
        a = linspace(0, 2*pi, 361);
        cx = CTRL.fcx.Value/1e3; cy = CTRL.fcy.Value/1e3;
        plot(ax, cx + CTRL.R.Value/1e3*sin(a), cy + CTRL.R.Value/1e3*cos(a), '--', 'Color', col);
        if CTRL.Rin.Value > 0
            plot(ax, cx + CTRL.Rin.Value/1e3*sin(a), cy + CTRL.Rin.Value/1e3*cos(a), ':', 'Color', col);
        end
        if abs(CTRL.azTo.Value - CTRL.azFrom.Value) < 360
            for az = [CTRL.azFrom.Value CTRL.azTo.Value]
                plot(ax, cx + [CTRL.Rin.Value CTRL.R.Value]/1e3*sind(az), cy + [CTRL.Rin.Value CTRL.R.Value]/1e3*cosd(az), '-', 'Color', col);
            end
        end
        plot(ax, P.ox/1e3, P.oy/1e3, 'p', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
    end

    function updateStatsLabel()
        if isempty(S.stats)
            lblStats.Text = 'No results yet.';
            return;
        end
        st = S.stats; dg = S.diag;
        if isempty(S.M), geo = '---';
        else
            Dm = HydrophoneGeometry.pairDist(S.M); Dn = Dm; Dn(1:size(Dm,1)+1:end) = inf;
            geo = sprintf('%d hydrophones, min spacing %.0f m, aperture %.0f m', size(S.M,1), min(Dn(:)), max(Dm(:)));
        end
        lblStats.Text = sprintf(['Pass %.1f%%  (score <= %.0f m, %s)\n' ...
            'Median %.1f m  |  Mean %.1f m  |  P95 %.1f m  |  Max %.1f m\n' ...
            'Detected %.1f%%  |  undetected %d, TDOA fail %d, solve fail %d  |  sources used %d / %d\n%s'], ...
            100*st.passRate, st.threshold_m, st.scoreMetric, st.medianErr, st.meanErr, st.p95Err, st.maxErr, ...
            100*st.detRate, dg.undetected, dg.tdoaFail, dg.solveFail, dg.gridUsed, dg.nP, geo);
    end

    function refreshAll()
        updateDepthLayers();
        refreshArrayTab();
        refreshEnvPlots();
        refreshGridPreview();
        refreshFieldPlot();
        updateStatsLabel();
    end

    function clearResults()
        S.stats = []; S.diag = struct();
        S.H = []; S.runInfo = [];
        S.sweep = [];
        if exist('btnReport', 'var') && isvalid(btnReport), btnReport.Enable = 'off'; end
        if exist('btnExportRun', 'var') && isvalid(btnExportRun), btnExportRun.Enable = 'off'; end
        if exist('lblStats', 'var') && isvalid(lblStats), updateStatsLabel(); end
    end

    function updateStatus(msg)
        lblStatus.Text = msg;
    end

%% ====================================================================
%  Evaluation, optimisation, sweep
%% ====================================================================
    function setBusy(b)
        S.busy = b;
        onoff = {'on','off'};
        st = onoff{1 + b};
        btnRun.Enable = st; btnOpt.Enable = st;
        if b, btnStop.Enable = 'on'; else, btnStop.Enable = 'off'; end
        if b, S.stopFlag = false; end
        drawnow;
    end

    function doStop()
        S.stopFlag = true;
        updateStatus('Stopping...');
    end

    function ok = confirmGridSize(n, what)
        ok = true;
        if n > 4000
            sel = uiconfirm(fig, sprintf('%s will simulate %d sources per evaluation. This may take a long time. Continue?', what, n), ...
                'Large grid', 'Options', {'Continue','Cancel'}, 'DefaultOption', 2, 'CancelOption', 2);
            ok = strcmp(sel, 'Continue');
        end
    end

    function doRunWAH()
        if S.busy, return; end
        if isempty(S.M), uialert(fig, 'No hydrophone geometry (Array Design tab).', 'Missing geometry'); return; end
        try
            P = readParams();
            ctx = buildContext(P);
            if ~confirmGridSize(size(ctx.grid,1), 'The evaluation'), return; end
            setBusy(true);
            S.progressPrefix = 'Evaluating:';
            t0 = tic;
            [st, dg] = evaluateGeometry(S.M, P.ox, P.oy, ctx);
            S.gridPts = ctx.grid; S.gridInfo = ctx.gridInfo;
            S.stats = st; S.diag = dg;
            initHistory(P, ctx);
            logIter(0, S.M, NaN, NaN, [], st, dg, scoreJ(st, P), true, 'WAH baseline');
            updateStatsLabel();
            if startsWith(CTRL.dispMetric.Value, 'Sweep'), CTRL.dispMetric.Value = 'Score error (m)'; end
            refreshFieldPlot();
            if dg.aborted
                updateStatus('Evaluation stopped - partial results shown.');
            else
                updateStatus(sprintf('Evaluation done in %.1f s: %.1f%% pass, median %.1f m.', toc(t0), 100*st.passRate, st.medianErr));
            end
            btnReport.Enable = 'off'; btnExportRun.Enable = 'on';
        catch ME
            uialert(fig, ME.message, 'Evaluation error');
            updateStatus('Evaluation failed');
        end
        setBusy(false);
    end

    function dirs = probeDirs(moveZ)
        ang = [0 180 90 270 45 225 135 315 22.5 202.5 112.5 292.5 67.5 247.5 157.5 337.5];
        dirs = [cosd(ang)' sind(ang)' zeros(numel(ang),1)];
        if moveZ
            dirs = [dirs(1:4,:); 0 0 1; 0 0 -1; dirs(5:end,:)];
        end
    end

    function doOptimiser()
        if S.busy, return; end
        if isempty(S.M), uialert(fig, 'No hydrophone geometry (Array Design tab).', 'Missing geometry'); return; end
        try
            P = readParams();
            ctx = buildContext(P);
            if ~confirmGridSize(size(ctx.grid,1), 'Each optimiser probe'), return; end
            [ok, msg] = HydrophoneGeometry.isValid(S.M, P.Lmax, P.minPair);
            if ~ok
                uialert(fig, sprintf('The starting geometry violates the constraints:\n%s\nAdjust Lmax / min spacing or the array.', msg), 'Invalid start');
                return;
            end
        catch ME
            uialert(fig, ME.message, 'Optimiser error');
            return;
        end

        setBusy(true);
        btnReport.Enable = 'off'; btnExportRun.Enable = 'off';
        S.H = []; S.runInfo = [];
        try
            S.progressPrefix = 'Baseline:';
            [bestStats, bestDiag] = evaluateGeometry(S.M, P.ox, P.oy, ctx);
            bestJ = scoreJ(bestStats, P);
            initHistory(P, ctx);
            logIter(0, S.M, NaN, NaN, [], bestStats, bestDiag, bestJ, true, 'baseline');
            acceptUpdate(bestStats, bestDiag, ctx);

            dirsAll = probeDirs(P.moveZ);
            N = size(S.M, 1);
            stall = 0;
            for it = 1:P.maxIters
                if S.stopFlag, break; end
                if P.maxIters > 1
                    stepFrac = max(P.minStep, P.initStep*(1 - (it-1)/(P.maxIters-1)));
                else
                    stepFrac = P.initStep;
                end
                stepSize = stepFrac * P.Lmax;
                hIdx = mod(it-1, N) + 1;
                dirs = dirsAll(1:min(P.probes, size(dirsAll,1)), :);
                cvec = mean(S.M(:,1:2), 1);
                rv = [S.M(hIdx,1:2) - cvec 0];
                if norm(rv) > 1e-9, rv = rv/norm(rv); else, rv = [1 0 0]; end
                dirs = [dirs; rv; -rv]; %#ok<AGROW>

                bestLocalJ = -inf; bestLocalM = []; bestLocalStats = []; bestLocalDiag = [];
                for k = 1:size(dirs,1)
                    if S.stopFlag, break; end
                    Mtry = S.M;
                    Mtry(hIdx,:) = Mtry(hIdx,:) + stepSize*dirs(k,:);
                    Mtry = HydrophoneGeometry.recentreXY(Mtry);
                    Mtry = HydrophoneGeometry.clampXY(Mtry, P.Lmax);
                    if strcmp(P.placeMode, 'seafloor'), Mtry(:,3) = max(Mtry(:,3), 0); end
                    if ~HydrophoneGeometry.isValid(Mtry, P.Lmax, P.minPair), continue; end
                    S.progressPrefix = sprintf('Iter %d/%d, H%d probe %d/%d:', it, P.maxIters, hIdx, k, size(dirs,1));
                    try
                        [stT, dgT] = evaluateGeometry(Mtry, P.ox, P.oy, ctx);
                    catch
                        continue;   % e.g. hydrophone moved onto land
                    end
                    if dgT.aborted, break; end
                    Jt = scoreJ(stT, P);
                    if Jt > bestLocalJ
                        bestLocalJ = Jt; bestLocalM = Mtry; bestLocalStats = stT; bestLocalDiag = dgT;
                    end
                end
                if S.stopFlag, break; end

                accepted = false;
                if ~isempty(bestLocalM) && bestLocalJ > bestJ + P.epsJ
                    S.M = bestLocalM; bestJ = bestLocalJ; bestStats = bestLocalStats; bestDiag = bestLocalDiag;
                    stall = 0; accepted = true; note = 'accepted (improved J)';
                    acceptUpdate(bestStats, bestDiag, ctx);
                    updateStatus(sprintf('Iter %d/%d: pass %.1f%%, median %.1f m, J = %.3f', it, P.maxIters, ...
                        100*bestStats.passRate, bestStats.medianErr, bestJ));
                else
                    stall = stall + 1;
                    note = sprintf('rejected (stall=%d)', stall);
                end
                logIter(it, S.M, hIdx, stepSize, dirs, bestStats, bestDiag, bestJ, accepted, note);

                if bestStats.passRate >= P.targetPass, break; end

                if stall >= P.stall
                    Msh = HydrophoneGeometry.shake(S.M, P.shake*P.Lmax, P.Lmax, P.minPair, P.moveZ);
                    if strcmp(P.placeMode, 'seafloor'), Msh(:,3) = max(Msh(:,3), 0); end
                    S.progressPrefix = 'Shake/reset:';
                    try
                        [stS, dgS] = evaluateGeometry(Msh, P.ox, P.oy, ctx);
                        S.M = Msh; bestStats = stS; bestDiag = dgS; bestJ = scoreJ(stS, P);
                        acceptUpdate(bestStats, bestDiag, ctx);
                        logIter(it + 0.5, S.M, NaN, NaN, [], bestStats, bestDiag, bestJ, true, 'shake/reset');
                    catch
                    end
                    stall = 0;
                end
                drawnow limitrate;
            end

            S.H.meta.timeEnd = datestr(now);
            S.H.meta.final = struct('M', S.M, 'stats', S.stats, 'diag', S.diag, 'J', scoreJ(S.stats, P));
            if S.stopFlag
                updateStatus(sprintf('Optimiser stopped: %.1f%% pass, median %.1f m.', 100*S.stats.passRate, S.stats.medianErr));
            else
                updateStatus(sprintf('Optimiser done: %.1f%% pass, median %.1f m.', 100*S.stats.passRate, S.stats.medianErr));
            end
        catch ME
            uialert(fig, ME.message, 'Optimiser error');
            updateStatus('Optimiser failed');
        end
        if ~isempty(S.H) && isfield(S.H, 'iter') && numel(S.H.iter) >= 1
            btnReport.Enable = 'on'; btnExportRun.Enable = 'on';
        end
        setBusy(false);
    end

    function acceptUpdate(st, dg, ctx)
        S.stats = st; S.diag = dg;
        S.gridPts = ctx.grid; S.gridInfo = ctx.gridInfo;
        refreshArrayTab();
        updateStatsLabel();
        if startsWith(CTRL.dispMetric.Value, 'Sweep'), CTRL.dispMetric.Value = 'Score error (m)'; end
        refreshFieldPlot();
    end

    function initHistory(P, ctx)
        S.H = struct();
        S.H.meta = struct();
        S.H.meta.version = 'WAH-i Ocean 1.0';
        S.H.meta.rngSeed = P.seed;
        S.H.meta.timeStart = datestr(now);
        S.H.meta.P = P;
        S.H.meta.env = ctx.env;
        S.H.meta.call = S.call;
        S.H.meta.band = S.band;
        S.H.meta.gridPts = ctx.grid;
        S.H.meta.gridInfo = ctx.gridInfo;
        S.H.meta.start = struct('source', 'GUI', 'M', S.M);
        S.H.iter = struct('it', {}, 'M', {}, 'hydIdx', {}, 'stepSize', {}, 'dirsTried', {}, ...
            'stats', {}, 'diag', {}, 'J', {}, 'accepted', {}, 'note', {});
        S.runInfo = struct('seed', P.seed, 'N', size(S.M,1), 'startNote', 'GUI start', 'mode', 'gui', ...
            'configTag', 'interactive', 'rep', 1, 'csvFile', '');
    end

    function logIter(it, M, hIdx, stepSize, dirs, st, dg, J, accepted, note)
        rec = struct('it', it, 'M', M, 'hydIdx', hIdx, 'stepSize', stepSize, 'dirsTried', dirs, ...
            'stats', st, 'diag', dg, 'J', J, 'accepted', accepted, 'note', char(note));
        S.H.iter(end+1) = rec;
    end

    function doSweep()
        if S.busy, return; end
        if isempty(S.M), uialert(fig, 'No hydrophone geometry.', 'Missing geometry'); return; end
        try
            P = readParams();
            ctx = buildContext(P);
            if strcmp(CTRL.sweepArea.Value, 'Around current origin')
                cx = P.ox; cy = P.oy; Rs = CTRL.sweepR.Value;
            else
                cx = P.fcx; cy = P.fcy; Rs = P.R;
            end
            stepS = CTRL.sweepStep.Value;
            n = floor(Rs/stepS);
            [Xo, Yo] = meshgrid((-n:n)*stepS, (-n:n)*stepS);
            keep = sqrt(Xo.^2 + Yo.^2) <= Rs + 1e-9;
            oxP = cx + Xo(keep); oyP = cy + Yo(keep);
            yawList = parseAngles(CTRL.sweepYaw.Value);
            if isempty(yawList), yawList = P.yaw; end
            nP = numel(oxP); nY = numel(yawList);
            ox = repmat(oxP, nY, 1); oy = repmat(oyP, nY, 1);
            yw = kron(yawList(:), ones(nP, 1));
            nO = numel(ox);
            if nO > 400
                error('%d candidates (%d positions x %d yaws) - increase the sweep spacing or use fewer yaws.', nO, nP, nY);
            end
            sel = uiconfirm(fig, sprintf('Evaluate %d candidates (%d positions x %d yaw angles) x %d sources = %d simulations?', ...
                nO, nP, nY, size(ctx.grid,1), nO*size(ctx.grid,1)), ...
                'Run sweep', 'Options', {'Run','Cancel'}, 'DefaultOption', 1, 'CancelOption', 2);
            if ~strcmp(sel, 'Run'), return; end
        catch ME
            uialert(fig, ME.message, 'Sweep error');
            return;
        end
        setBusy(true);
        sw = struct('ox', ox, 'oy', oy, 'yaw', yw, 'pass', nan(nO,1), 'median', nan(nO,1), 'mean', nan(nO,1), ...
            'p95', nan(nO,1), 'det', nan(nO,1), 'valid', false(nO,1), 'best', NaN, 'rank', CTRL.sweepRank.Value);
        try
            for q = 1:nO
                if S.stopFlag, break; end
                S.progressPrefix = sprintf('Sweep candidate %d/%d (yaw %.0f):', q, nO, yw(q));
                try
                    [stq, dgq] = evaluateGeometry(S.M, ox(q), oy(q), ctx, yw(q));
                    if dgq.aborted, break; end
                    sw.pass(q) = stq.passRate; sw.median(q) = stq.medianErr; sw.mean(q) = stq.meanErr;
                    sw.p95(q) = stq.p95Err; sw.det(q) = stq.detRate; sw.valid(q) = true;
                catch
                    sw.valid(q) = false;  % on land etc.
                end
            end
            if any(sw.valid)
                idx = find(sw.valid);
                if strcmp(CTRL.sweepRank.Value, 'Pass rate')
                    key = [-sw.pass(idx) sw.median(idx)];
                else
                    key = [sw.median(idx) -sw.pass(idx)];
                end
                [~, ord] = sortrows(key);
                sw.best = idx(ord(1));
                lblSweep.Text = sprintf('Best of %d candidates: origin (%.0f, %.0f) m, yaw %.1f deg - pass %.1f%%, median %.1f m.', ...
                    sum(sw.valid), sw.ox(sw.best), sw.oy(sw.best), sw.yaw(sw.best), 100*sw.pass(sw.best), sw.median(sw.best));
            else
                lblSweep.Text = 'No valid positions (all on land?).';
            end
            S.sweep = sw;
            CTRL.dispMetric.Value = 'Sweep: pass rate (%)';
            refreshFieldPlot();
            updateStatus(sprintf('Sweep finished (%d candidates evaluated).', sum(sw.valid)));
        catch ME
            uialert(fig, ME.message, 'Sweep error');
        end
        setBusy(false);
    end

    function doApplyBestSweep()
        if isempty(S.sweep) || ~isfinite(S.sweep.best), updateStatus('Run a sweep first.'); return; end
        CTRL.ox.Value = S.sweep.ox(S.sweep.best);
        CTRL.oy.Value = S.sweep.oy(S.sweep.best);
        CTRL.yaw.Value = S.sweep.yaw(S.sweep.best);
        sw = S.sweep;
        clearResults();
        S.sweep = sw;   % keep the sweep visible
        refreshArrayTab(); refreshEnvPlots(); refreshFieldPlot();
        updateStatus(sprintf('Array moved to (%.0f, %.0f) m, yaw %.1f deg. Evaluate or optimise there.', CTRL.ox.Value, CTRL.oy.Value, CTRL.yaw.Value));
    end

    function doExportSweep()
        if isempty(S.sweep), updateStatus('Run a sweep first.'); return; end
        [file, pth] = uiputfile({'*.csv','CSV (*.csv)'}, 'Export sweep', 'array_position_sweep.csv');
        figure(fig);
        if isequal(file, 0), return; end
        try
            sw = S.sweep;
            fid = fopen(fullfile(pth, file), 'w');
            fprintf(fid, 'origin_x_m,origin_y_m,yaw_deg,valid,pass_rate,median_err_m,mean_err_m,p95_err_m,detected_frac\n');
            for q = 1:numel(sw.ox)
                fprintf(fid, '%.2f,%.2f,%.3f,%d,%.4f,%.3f,%.3f,%.3f,%.4f\n', sw.ox(q), sw.oy(q), sw.yaw(q), sw.valid(q), ...
                    sw.pass(q), sw.median(q), sw.mean(q), sw.p95(q), sw.det(q));
            end
            fclose(fid);
            updateStatus(sprintf('Sweep saved to %s', file));
        catch ME
            uialert(fig, ME.message, 'Export error');
        end
    end

%% ====================================================================
%  Inspect a single source
%% ====================================================================
    function doInspect()
        if isempty(S.M), uialert(fig, 'No hydrophone geometry.', 'Missing geometry'); return; end
        d = uifigure('Name', 'Inspect a source', 'Position', centreOnParent([900 620]));
        g = uigridlayout(d, [3 1]);
        g.RowHeight = {30, '1x', 60};
        gtop = uigridlayout(g, [1 8]);
        gtop.Padding = [0 0 0 0];
        gtop.ColumnWidth = {'fit', 90, 'fit', 90, 'fit', 70, 110, '1x'};
        uilabel(gtop, 'Text', 'East X (m):');
        exX = uieditfield(gtop, 'numeric', 'Value', CTRL.fcx.Value + CTRL.R.Value/2);
        uilabel(gtop, 'Text', 'North Y (m):');
        exY = uieditfield(gtop, 'numeric', 'Value', CTRL.fcy.Value);
        uilabel(gtop, 'Text', 'Depth (m):');
        dList = parseDepths(CTRL.depths.Value);
        exD = uieditfield(gtop, 'numeric', 'Value', dList(1));
        uibutton(gtop, 'Text', 'Simulate', 'ButtonPushedFcn', @(~,~)runInspect());
        axI = uiaxes(g);
        lbI = uilabel(g, 'Text', '', 'WordWrap', 'on', 'VerticalAlignment', 'top');
        runInspect();

        function runInspect()
            try
                P = readParams();
                ctx = buildContext(P);
                Hn = HydrophoneGeometry.toAbsolute(S.M, placement(P), ctx.env);
                src = [exX.Value exY.Value -abs(exD.Value)];
                cfg = wahi_geometry_cfg(ctx.cfg, P, ctx.env, Hn, src);
                loc = WhaleCallLocaliser(cfg, ctx.env, Hn);
                rng(P.seed, 'twister');
                sim = loc.simulate(src);
                cla(axI); hold(axI, 'on'); grid(axI, 'on');
                [snrAll, rAll] = loc.receivedSNR(src, loc.Htrue);
                if isempty(sim.signals)
                    title(axI, sprintf('Not detected: only %d hydrophone(s) above DT', numel(sim.idx)));
                    lbI.Text = sprintf('Received SNR per hydrophone: %s dB (DT = %.1f dB). Ranges: %s km', ...
                        mat2str(round(max(snrAll,-99),1)), P.DT, mat2str(round(rAll/1e3,2)));
                    hold(axI, 'off');
                    return;
                end
                sig = sim.signals;
                tt = (0:size(sig,1)-1)/sim.fs;
                sc = 0.45/max(abs(sig(:)) + eps);
                for q = 1:size(sig,2)
                    plot(axI, tt, q + sc*sig(:,q), 'Color', [0.1 0.3 0.6]);
                end
                yticks(axI, 1:size(sig,2));
                yticklabels(axI, compose('H%d', sim.idx));
                xlabel(axI, 'Time relative to first arrival (s)');
                err = norm(sim.est(1:2) - src(1:2));
                title(axI, sprintf('Received signals (band-limited, with noise). Horizontal error = %.1f m', err));
                hold(axI, 'off');
                lbI.Text = sprintf(['True TDOA (s): %s\nEst. TDOA (s): %s\nSNR (dB): %s | estimate (%.1f, %.1f, %.1f) m ' ...
                    '| true (%.1f, %.1f, %.1f) m'], mat2str(sim.trueTdoa, 5), mat2str(sim.tdoa, 5), ...
                    mat2str(round(sim.snr(sim.idx), 1)), sim.est(1), sim.est(2), sim.est(3), src(1), src(2), src(3));
            catch ME
                lbI.Text = ['Error: ' ME.message];
            end
        end
    end

    function pos = centreOnParent(sz)
        p = fig.Position;
        pos = [round(p(1) + (p(3)-sz(1))/2), round(p(2) + (p(4)-sz(2))/2), sz(1), sz(2)];
    end

%% ====================================================================
%  Exports, report, sessions
%% ====================================================================
    function doExportResults()
        if isempty(S.stats), uialert(fig, 'Evaluate the field first.', 'No results'); return; end
        [file, pth] = uiputfile({'*.csv','CSV (*.csv)'}, 'Export per-source results', ...
            sprintf('wahi_results_%s.csv', datestr(now, 'yyyymmdd_HHMMSS')));
        figure(fig);
        if isequal(file, 0), return; end
        try
            st = S.stats; G = S.gridPts;
            fid = fopen(fullfile(pth, file), 'w');
            fprintf(fid, 'x_m,y_m,depth_m,flag,used,pass,score_err_m,horiz_err_m,err3d_m,depth_err_m,bearing_err_deg,n_detecting,min_snr_db,max_snr_db,est_x_m,est_y_m,est_depth_m\n');
            for q = 1:size(G,1)
                fprintf(fid, '%.2f,%.2f,%.2f,%d,%d,%d,%.3f,%.3f,%.3f,%.3f,%.4f,%d,%.2f,%.2f,%.2f,%.2f,%.2f\n', ...
                    G(q,1), G(q,2), -G(q,3), st.flag(q), st.usedMask(q), st.passMask(q), st.score_m(q), st.errH_m(q), ...
                    st.err3_m(q), st.errZ_m(q), st.brg_deg(q), st.nDet(q), st.minSNR(q), st.maxSNR(q), ...
                    st.est(q,1), st.est(q,2), -st.est(q,3));
            end
            fclose(fid);
            updateStatus(sprintf('Results saved to %s (flag: 1 ok, 2 undetected, 3 TDOA fail, 4 solve fail, 0 skipped)', file));
        catch ME
            uialert(fig, ME.message, 'Export error');
        end
    end

    function doExportRun()
        if isempty(S.H), uialert(fig, 'No run history. Evaluate or optimise first.', 'Nothing to export'); return; end
        [file, pth] = uiputfile({'*.mat','MAT-files (*.mat)'}, 'Save run as .mat', ...
            sprintf('wahi_ocean_run__%dhyd__%s.mat', size(S.M,1), datestr(now, 'yyyymmdd_HHMMSS')));
        figure(fig);
        if isequal(file, 0), updateStatus('Export cancelled'); return; end
        try
            H = S.H; runInfo = S.runInfo; %#ok<NASGU>
            save(fullfile(pth, file), 'H', 'runInfo', '-v7.3');
            updateStatus(sprintf('Saved run: %s', file));
        catch ME
            uialert(fig, ME.message, 'Export error');
        end
    end

    function doExportWorkspace()
        if isdeployed, return; end
        assignin('base', 'wahi_M', S.M);
        assignin('base', 'wahi_gridPts', S.gridPts);
        assignin('base', 'wahi_stats', S.stats);
        assignin('base', 'wahi_diag', S.diag);
        assignin('base', 'wahi_H', S.H);
        try, assignin('base', 'wahi_env', currentEnv(readParams(false))); catch, assignin('base', 'wahi_env', S.env); end
        assignin('base', 'wahi_call', S.call);
        assignin('base', 'wahi_sweep', S.sweep);
        try, assignin('base', 'wahi_P', readParams()); catch, end
        updateStatus('Exported to workspace: wahi_M, wahi_gridPts, wahi_stats, wahi_diag, wahi_H, wahi_env, wahi_call, wahi_sweep, wahi_P');
    end

    function doReport()
        if isempty(S.H) || numel(S.H.iter) < 1
            uialert(fig, 'No run history. Optimise first.', 'No data'); return;
        end
        try
            f = wahi_ocean_report(S.H);
            movegui(f, 'center');
            choice = uiconfirm(fig, 'Save report figure as PNG?', 'Save report?', ...
                'Options', {'Yes','No'}, 'DefaultOption', 2, 'CancelOption', 2);
            if strcmp(choice, 'Yes')
                [file, pth] = uiputfile({'*.png','PNG (*.png)'}, 'Save report', 'wahi_ocean_report.png');
                if ~isequal(file, 0)
                    exportgraphics(f, fullfile(pth, file), 'Resolution', 180);
                    updateStatus(sprintf('Report saved: %s', file));
                end
            else
                updateStatus('Report generated');
            end
        catch ME
            uialert(fig, ME.message, 'Report error');
        end
    end

    function doSaveSession()
        [file, pth] = uiputfile({'*.mat','WAH-i Ocean session (*.mat)'}, 'Save session', ...
            sprintf('wahi_ocean_session_%s.mat', datestr(now, 'yyyymmdd_HHMMSS')));
        figure(fig);
        if isequal(file, 0), return; end
        try
            session = struct();
            session.version = 'WAH-i Ocean 1.0';
            keys = fieldnames(CTRL);
            vals = struct();
            for q = 1:numel(keys), vals.(keys{q}) = CTRL.(keys{q}).Value; end
            session.values = vals;
            session.M = S.M;
            session.bathy = S.env.bathy;
            session.lat0 = S.env.lat0; session.lon0 = S.env.lon0;
            session.sspProfile = S.sspProfile;
            session.ncom = S.ncom;
            session.callTemplate = S.callTemplate;
            session.tplLabel = lblTpl.Text;
            session.sspLabel = lblSSP.Text; %#ok<STRNU>
            save(fullfile(pth, file), 'session', '-v7.3');
            updateStatus(sprintf('Session saved: %s', file));
        catch ME
            uialert(fig, ME.message, 'Save error');
        end
    end

    function doLoadSession()
        [file, pth] = uigetfile({'*.mat','WAH-i Ocean session (*.mat)'}, 'Load session');
        figure(fig);
        if isequal(file, 0), return; end
        try
            L = load(fullfile(pth, file));
            if ~isfield(L, 'session'), error('Not a WAH-i Ocean session file.'); end
            ss = L.session;
            keys = fieldnames(ss.values);
            for q = 1:numel(keys)
                if isfield(CTRL, keys{q})
                    try
                        if strcmp(keys{q}, 'dispLayer'), continue; end
                        val = ss.values.(keys{q});
                        if strcmp(keys{q}, 'placeMode'), val = placeLabel(placeCode(val)); end
                        CTRL.(keys{q}).Value = val;
                    catch
                    end
                end
            end
            S.M = ss.M;
            S.env.bathy = ss.bathy;
            if ~isempty(S.env.bathy), S.env.bathy = OceanEnvironment.finaliseBathy(S.env.bathy); end
            S.env.lat0 = ss.lat0; S.env.lon0 = ss.lon0;
            S.sspProfile = ss.sspProfile;
            S.ncom = ss.ncom;
            if ~isempty(S.ncom), S.ncom = OceanEnvironment.rebuildField(S.ncom); end
            S.callTemplate = ss.callTemplate;
            if isfield(ss, 'tplLabel'), lblTpl.Text = ss.tplLabel; end
            if isfield(ss, 'sspLabel'), lblSSP.Text = ss.sspLabel; end
            setPresetNote();
            updateBathyInfo();
            updateCall(true);
            clearResults(); refreshAll();
            updateStatus(sprintf('Session loaded: %s', file));
        catch ME
            uialert(fig, ME.message, 'Load error');
        end
    end

    function doResetDefaults()
        choice = uiconfirm(fig, 'Reset all settings to defaults? (Loaded bathymetry, profiles and templates are kept.)', ...
            'Reset defaults', 'Options', {'Reset','Cancel'}, 'DefaultOption', 2, 'CancelOption', 2);
        if ~strcmp(choice, 'Reset'), return; end
        keys = fieldnames(DEF);
        for q = 1:numel(keys)
            try
                if strcmp(keys{q}, 'dispLayer'), continue; end
                CTRL.(keys{q}).Value = DEF.(keys{q});
            catch
            end
        end
        applyPreset(true);
        S.M = HydrophoneGeometry.pattern(CTRL.patName.Value, patternOpts());
        onOriginChanged();
        clearResults(); refreshAll();
        updateStatus('Defaults restored.');
    end

    function showTips()
        txt = sprintf([ ...
            'Workflow\n' ...
            '1) Array Design: generate a pattern (dx/dy spacing), then edit cells for uneven layouts. Export CSV to keep it.\n' ...
            '2) Environment: optionally load bathymetry (GEBCO NetCDF, CSV x/y/z or lon/lat/z, MAT) and a sound-speed\n' ...
            '   profile or NCOM file. Check the SNR-vs-range plot: sources beyond the detection range will be undetected.\n' ...
            '3) Calls: pick a species preset or import a template. Keep the sample rate as low as the call allows -\n' ...
            '   long-range runs scale with call length x sample rate.\n' ...
            '4) Field & Optimise: evaluate, then optimise. Use "Array position sweep" to find a good deployment site.\n\n' ...
            'Speed\n' ...
            '- Each source costs a few ms; the optimiser evaluates the whole field (probes+2) times per iteration.\n' ...
            '  Optimise on a coarse spacing first, then evaluate the result on a fine grid.\n' ...
            '- Enable the parallel pool for large grids.\n\n' ...
            'Interpreting results\n' ...
            '- Range error grows quickly beyond a few apertures; bearing error stays small. Look at "Bearing error".\n' ...
            '- Set the localiser sound speed different from the truth model to see the effect of SSP mismatch.\n' ...
            '- Set the assumed source depth different from the true depth to see depth-mismatch errors in 2-D solves.\n' ...
            '- Clock error of 1 ms ~ 1.5 m of range difference; at long range this is amplified strongly.\n\n' ...
            'Refraction / ducting\n' ...
            '- Environment > Propagation model = Ray tracing: Snell''s-law rays through the profile at the array,\n' ...
            '  with surface (-1) and bottom (BL dB) reflections. SNR and arrivals then show shadow zones,\n' ...
            '  convergence zones and multipath. Set "Localiser assumes" = Ray-traced so the solver also\n' ...
            '  accounts for refraction (otherwise you see the error of a straight-ray localiser).\n\n' ...
            'Array orientation\n' ...
            '- Yaw / pitch / roll rotate the designed frame about the array origin (roll -> pitch -> yaw).\n' ...
            '  The World coordinates table shows where each hydrophone ends up (E, N, height, depth).\n' ...
            '- The sweep can try several yaw angles at each position.\n\n' ...
            'Limitations\n' ...
            '- Ray model is range-independent (one profile and a flat bottom at the array), geometric acoustics\n' ...
            '  only (no diffraction into shadow zones, no sediment model; incoherent TL).\n' ...
            '- No source directivity.\n' ...
            '- Species presets are simple parametric approximations - check against your recordings.\n']);
        d = uifigure('Name', 'Tips', 'Position', centreOnParent([780 640]));
        g = uigridlayout(d, [2 1]); g.RowHeight = {'1x', 34};
        uitextarea(g, 'Value', strsplit(txt, newline), 'Editable', 'off');
        uibutton(g, 'Text', 'Close', 'ButtonPushedFcn', @(~,~)delete(d));
    end

    function showAbout()
        txt = sprintf([ ...
            'WAH-i Ocean - hydrophone array design & optimisation for whale localisation\n\n' ...
            'Modified from WAH-i (Widefield Acoustics Heuristic - inverse, iterative)\n' ...
            'Original author: Ravi Umadi (2026), https://github.com/raviumadi/WAHi-Optimiser\n' ...
            'Original paper: Umadi (2026) WAH-i: Optimising Microphone Array Geometry for Customised\n' ...
            'Localisation Accuracy. https://doi.org/10.64898/2026.02.07.704547\n\n' ...
            'This version replaces the in-air bat-call model with an underwater model: hydrophones,\n' ...
            'ray tracing (refraction, ducting, multipath), 3-D array orientation,\n' ...
            'bathymetry, sound-speed profiles (constant / 1-D / NCOM 3-D), sonar-equation SNR,\n' ...
            'detection thresholds, clock and survey uncertainty, whale call presets and templates,\n' ...
            'fields of accuracy to 100 km and an array-position sweep.\n\n' ...
            'Licence: GPLv3. Provided as-is without warranty; validate results for your deployment.\n']);
        d = uifigure('Name', 'About', 'Position', centreOnParent([700 380]));
        g = uigridlayout(d, [2 1]); g.RowHeight = {'1x', 34};
        uitextarea(g, 'Value', strsplit(txt, newline), 'Editable', 'off');
        uibutton(g, 'Text', 'Close', 'ButtonPushedFcn', @(~,~)delete(d));
    end
end

%% ===================== Local (non-nested) helpers =====================

function resetAxesModes(ax)
ax.XLimMode = 'auto'; ax.YLimMode = 'auto'; ax.ZLimMode = 'auto';
ax.DataAspectRatioMode = 'auto'; ax.PlotBoxAspectRatioMode = 'auto';
ax.CLimMode = 'auto';
end

function [Xc, Yc, Dc] = bathyWindow(B, xr, yr)
% Bathymetry sub-grid covering x/y ranges (world metres) for contouring
nx = 80; ny = 80;
xv = linspace(xr(1), xr(2), nx); yv = linspace(yr(1), yr(2), ny);
[Xc, Yc] = meshgrid(xv, yv);
Dc = B.F(Yc, Xc);
end

function [xd, yd, Dd] = decimateForPlot(B, nmax)
sx = max(1, ceil(numel(B.x)/nmax)); sy = max(1, ceil(numel(B.y)/nmax));
xd = B.x(1:sx:end); yd = B.y(1:sy:end);
Dd = B.depth(1:sy:end, 1:sx:end);
end
