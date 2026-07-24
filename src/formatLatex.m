function formatLatex(ax)
% Applies LaTeX formatting and consistent styling to a given axes handle.
% Usage: formatAxesLatex(gca)

    set(ax, 'TickLabelInterpreter', 'latex', ...
            'FontSize', 10, ...
            'FontWeight', 'bold');
    
    grid(ax, 'on');
    grid(ax, 'minor');
    
    % If labels or title are already set, preserve and format them
    if ~isempty(get(get(ax, 'XLabel'), 'String'))
        set(get(ax, 'XLabel'), 'Interpreter', 'latex', ...
                               'FontSize', 10, ...
                               'FontWeight', 'bold');
    end

    if ~isempty(get(get(ax, 'YLabel'), 'String'))
        set(get(ax, 'YLabel'), 'Interpreter', 'latex', ...
                               'FontSize', 10, ...
                               'FontWeight', 'bold');
    end

    if ~isempty(get(get(ax, 'Title'), 'String'))
        set(get(ax, 'Title'), 'Interpreter', 'latex', ...
                              'FontSize', 12, ...
                              'FontWeight', 'bold');
    end
    
    if ~isempty(get(get(ax, 'Subtitle'), 'String'))
        set(get(ax, 'Subtitle'), 'Interpreter', 'latex', ...
                              'FontSize', 10, ...
                              'FontWeight', 'bold');
    end
end