function imx662_table_fig(name, ttl, headers, data, ~)
%IMX662_TABLE_FIG  Draw a results table as a figure in booktabs style:
%   heavy top rule, thin rule under the header, heavy bottom rule, no
%   vertical lines, no shading. (Publication versions of every table are
%   also built in LaTeX: see 11_Summary_tables/*.tex.)
    [m, n] = size(data);
    f = figure('Name', name, 'Color', 'w', 'Position', [40 40 min(420 + 170*n, 2000) 150 + 34*m]);
    ax = axes(f, 'Position', [0.02 0.04 0.96 0.80]); axis(ax, 'off'); hold(ax, 'on');
    xlim(ax, [0 n]); ylim(ax, [0 m + 1.3]);
    w = ones(1, n); w(1) = 2.4; x = [0 cumsum(w)] / sum(w) * n;
    fn = 'Times New Roman';
    for j = 1:n
        if j == 1, xa = x(j) + 0.05; ha = 'left'; else, xa = x(j+1) - 0.05; ha = 'right'; end
        text(ax, xa, m + 0.62, headers{j}, 'HorizontalAlignment', ha, 'FontSize', 11, ...
             'FontName', fn, 'Interpreter', 'none');
    end
    for i = 1:m
        for j = 1:n
            v = data{i, j};
            if isnumeric(v)
                if isempty(v) || isnan(v), s = '---'; elseif abs(v) >= 1e5 || (abs(v) < 1e-3 && v ~= 0), s = sprintf('%.3g', v);
                elseif v == round(v), s = sprintf('%d', round(v)); else, s = sprintf('%.4g', v); end
            else
                s = char(strrep(string(v), 'This analysis', 'Carissma (2026)'));
            end
            if j == 1, xa = x(j) + 0.05; ha = 'left'; else, xa = x(j+1) - 0.05; ha = 'right'; end
            text(ax, xa, m - i + 0.5, s, 'HorizontalAlignment', ha, 'FontSize', 11, 'FontName', fn, 'Interpreter', 'none');
        end
    end
    plot(ax, [0 n], [m+1.25 m+1.25], 'k-', 'LineWidth', 1.4);   % \toprule
    plot(ax, [0 n], [m m], 'k-', 'LineWidth', 0.6);             % \midrule
    plot(ax, [0 n], [0 0], 'k-', 'LineWidth', 1.4);             % \bottomrule
    title(ax, ttl, 'FontSize', 12, 'FontName', fn, 'FontWeight', 'normal', 'Interpreter', 'none', ...
          'HorizontalAlignment', 'left', 'Position', [0 m + 1.35 0]);
end
