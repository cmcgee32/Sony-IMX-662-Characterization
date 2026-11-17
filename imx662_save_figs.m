function imx662_save_figs(prefix, outDir)
% Save every open figure as PNG (for reports) and .fig (editable in MATLAB), then close.
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    figs = findobj(0, 'Type', 'figure');
    [~, o] = sort([figs.Number]); figs = figs(o);
    for i = 1:numel(figs)
        nm = regexprep(figs(i).Name, '[^A-Za-z0-9]+', '_');
        base = fullfile(outDir, sprintf('%s_%02d_%s', prefix, i, nm));
        try
            exportgraphics(figs(i), [base '.png'], 'Resolution', 130);
        catch ME
            warning('PNG export failed for %s: %s', base, ME.message);
        end
        try, savefig(figs(i), [base '.fig']); catch, end
    end
    close all;
end
