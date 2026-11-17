function imx662_style(fig)
%IMX662_STYLE  Clean, publication-style look for every axes in a figure
%   (in the spirit of the gramm toolbox: white panels, light grey grid,
%   outward ticks, no box, unboxed legends, consistent fonts).
    if nargin < 1, fig = gcf; end
    set(fig, 'Color', 'w');
    ax = findobj(fig, 'Type', 'axes');
    for a = ax'
        if strcmp(a.Visible, 'off'), continue; end
        set(a, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.008 0.008], 'LineWidth', 0.9, ...
            'FontName', 'Helvetica', 'FontSize', 11, 'XColor', [0.15 0.15 0.15], ...
            'GridColor', [0.80 0.80 0.80], 'GridAlpha', 0.7, 'MinorGridAlpha', 0.25, 'Layer', 'top', 'Color', 'w');
        if numel(a.YAxis) == 1, a.YColor = [0.15 0.15 0.15]; end
        grid(a, 'on');
        a.XMinorGrid = 'off'; a.YMinorGrid = 'off';
        if ~isempty(findobj(a, 'Type', 'image')), grid(a, 'off'); end   % no grid lines over images
        a.Title.FontWeight = 'bold'; a.Title.FontSize = 13;
        a.XLabel.FontSize = 12; a.YLabel.FontSize = 12;
        % right y-axis (yyaxis) keeps its own colour
    end
    lg = findobj(fig, 'Type', 'legend');
    for l = lg', set(l, 'Box', 'off', 'FontSize', 10, 'Color', 'none'); end
    cb = findobj(fig, 'Type', 'colorbar');
    for b = cb', set(b, 'TickDirection', 'out', 'Box', 'off', 'FontSize', 10); end
end
