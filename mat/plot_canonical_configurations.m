%% plot_canonical_configurations.m
% Plot the six canonical microphone-array starting configurations used in
% the WAH-i optimisation study, labelled with configuration IDs 1--6.
%
% The script reads the authoritative coordinate CSV files in configs/
% and exports manuscript-ready PDF and PNG figures to fig/.

clear; clc;

%% User options
figureSizePixels = [750 450];  % [width height]
saveFigure = true;            % true exports PDF and PNG

%% Paths
scriptDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(scriptDir);
configDir = fullfile(projectDir, 'configs');
figureDir = fullfile(projectDir, 'fig');

if ~isfolder(figureDir)
    mkdir(figureDir);
end

%% Authoritative configuration-ID mapping (Figure 1 and Table 2)
configs = struct( ...
    'id',       num2cell(1:6), ...
    'name',     {'Pyramid', 'Tetrahedron', 'Octahedron', ...
                 'Double tetrahedron', 'Pentagonal antiprism', 'Icosahedron'}, ...
    'file',     {'4mics_Pyramid.csv', ...
                 '4mics_Tetrahedron.csv', ...
                 '6mics_Octahedron.csv', ...
                 '8mics_double_tetrahedorn.csv', ...
                 '10mics_pentagonal_antiprism.csv', ...
                 '12mics_icosahedron.csv'} ...
    );

%% Load and validate coordinates
allCoordinates = [];

for k = 1:numel(configs)
    csvPath = fullfile(configDir, configs(k).file);
    if ~isfile(csvPath)
        error('Missing configuration file: %s', csvPath);
    end

    coordinates = readmatrix(csvPath, 'FileType', 'text');
    coordinates = coordinates(:, 1:min(3, size(coordinates, 2)));
    coordinates = coordinates(all(isfinite(coordinates), 2), :);

    if size(coordinates, 2) ~= 3 || size(coordinates, 1) < 4
        error('Expected at least four rows of x,y,z coordinates in: %s', csvPath);
    end

    configs(k).coordinates = coordinates;
    configs(k).nMicrophones = size(coordinates, 1);
    allCoordinates = [allCoordinates; coordinates]; %#ok<AGROW>
end

% Use identical axis limits in every panel so array sizes remain comparable.
axisLimit = 1.12 * max(abs(allCoordinates), [], 'all');

%% Figure styling
colours.microphone = [0.10 0.38 0.72];
colours.edge = [0.25 0.28 0.32];
colours.label = [0.08 0.08 0.08];

fig = figure( ...
    'Name', 'Canonical microphone-array configurations', ...
    'Color', 'white', ...
    'Units', 'pixels', ...
    'Position', [100 100 figureSizePixels]);

layout = tiledlayout(fig, 2, 3, ...
    'Padding', 'compact', ...
    'TileSpacing', 'compact');

for k = 1:numel(configs)
    ax = nexttile(layout, k);
    hold(ax, 'on');

    M = configs(k).coordinates;
    edges = convexHullEdges(M);

    for e = 1:size(edges, 1)
        pair = edges(e, :);
        plot3(ax, M(pair, 1), M(pair, 2), M(pair, 3), '-', ...
            'Color', colours.edge, ...
            'LineWidth', 1.0);
    end

    scatter3(ax, M(:, 1), M(:, 2), M(:, 3), 24, ...
        'MarkerFaceColor', colours.microphone, ...
        'MarkerEdgeColor', 'white', ...
        'LineWidth', 0.8);

    % Number microphones to make the coordinate layouts unambiguous.
    labelOffset = 0.035 * axisLimit;
    for microphone = 1:size(M, 1)
        text(ax, M(microphone, 1) + labelOffset, ...
            M(microphone, 2) + labelOffset, ...
            M(microphone, 3) + labelOffset, ...
            string(microphone), ...
            'Color', colours.label, ...
            'FontSize', 7, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'bottom');
    end

    title(ax, sprintf('ID %d: %s (N = %d)', ...
        configs(k).id, configs(k).name, configs(k).nMicrophones), ...
        'FontWeight', 'bold', ...
        'FontSize', 9);
    
    xlabel(ax, 'x (m)');
    ylabel(ax, 'y (m)');
    zlabel(ax, 'z (m)', 'Interpreter', 'latex');
    xlim(ax, [-axisLimit axisLimit]);
    ylim(ax, [-axisLimit axisLimit]);
    zlim(ax, [-axisLimit axisLimit]);
    daspect(ax, [1 1 1]);
    view(ax, 38, 24);
    % grid(ax, 'on');
    box(ax, 'on');
    formatLatex(gca);
    ax.GridAlpha = 0.12;
    ax.MinorGridAlpha = 0.06;
    ax.FontSize = 7;
    ax.LineWidth = 0.7;
end

title(layout, 'Canonical starting configurations used for optimisation', ...
    'FontSize', 12, ...
    'FontWeight', 'bold', 'Interpreter', 'latex');

%% Export
if saveFigure
    pdfPath = fullfile(figureDir, 'canonical_configurations.pdf');
    pngPath = fullfile(figureDir, 'canonical_configurations.png');

    exportgraphics(fig, pdfPath, ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'white');
    exportgraphics(fig, pngPath, ...
        'Resolution', 300, ...
        'BackgroundColor', 'white');

    fprintf('Saved canonical-configuration figures:\n  %s\n  %s\n', ...
        pdfPath, pngPath);
end

%% Local functions
function edges = convexHullEdges(coordinates)
% Return the unique line segments forming the 3-D convex hull.
    hullFaces = convhulln(coordinates);
    faceEdges = [hullFaces(:, [1 2]); ...
                 hullFaces(:, [2 3]); ...
                 hullFaces(:, [3 1])];
    edges = unique(sort(faceEdges, 2), 'rows');
end
