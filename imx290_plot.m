function pg = imx290_plot(E, OUT)
%IMX290_PLOT  Figures for the IMX290 per-pixel gain calibration (from imx290_ecc_analysis).
%   imx290_plot(E)        E = struct returned/saved by imx290_ecc_analysis
%   imx290_plot()         reloads Downloads/IMX290_results/imx290_results.mat and redraws
    c = imx662_cfg();
    if nargin < 2, OUT = fullfile(c.D, 'IMX290_results'); end
    if nargin < 1 || isempty(E), E = load(fullfile(OUT, 'imx290_results.mat')); end
    if ~isfield(E, 'nAll'), E.nAll = 1691; end
    g = E.g; G = E.G; Fr = E.Fr; Fc = E.Fc; Fr2 = E.Fr2; Fc2 = E.Fc2; acc = E.acc; ben = E.ben; XY = E.XY;
    fanoMn = 2.3548*sqrt(c.fano*c.epsSi*c.E.MnKa);
    gm = @(p, x) p(1)*exp(-0.5*((x-p(2))/p(3)).^2) + p(4);
    gpk = @(p, x) p(1)*exp(-0.5*((x-p(2))/p(3)).^2);
    % spectra for display: 16 eV bins (wider than the native 13.6 eV spacing, so no comb pattern)
    w = 16; eb = 0:w:8000; e16 = eb(1:end-1) + w/2;
    c16 = accumarray(max(1, ceil(E.ec(:)/w)), E.cor(:), [numel(e16) 1])';
    r16 = accumarray(max(1, min(numel(e16), ceil(E.erw(:)/w))), E.raw(:), [numel(e16) 1])';
    sR = w/4;   % scale of the 4-eV-bin fits to 16-eV bins (per-pixel spectrum); raw fit is per native bin
    sRaw = w/(4/G);

    % 1. spectrum before / after
    f = figure('Name', 'IMX290 Fe-55 spectrum before and after gain correction', 'Position', [40 40 1300 700]);
    stairs(e16, max(r16, 0.5), 'Color', [0.6 0.6 0.6], 'LineWidth', 1.3); hold on;
    stairs(e16, max(c16, 0.5), 'Color', c.col.n1, 'LineWidth', 1.3);
    for q = {Fc, Fc2}
        p = q{1}.p; xx = linspace(p(2) - 2.5*p(3), p(2) + 2.5*p(3), 300);
        plot(xx, (gpk(p, xx) + max(p(4), 0))*sR, '--', 'Color', c.col.fit, 'LineWidth', 1.5);
    end
    set(gca, 'YScale', 'log'); xlim([1000 7500]); ylim([100 max(c16)*20]);
    legend({'Single gain (median)', 'Per-pixel ECC gain', 'Gaussian Fit'}, 'Location', 'northwest');
    text(0.62, 0.95, sprintf(['Mn K\\alpha FWHM\nSingle gain: %.1f \\pm %.1f eV\nPer-pixel gain: %.1f \\pm %.1f eV\n\n' ...
         'Mn K\\beta FWHM\nSingle gain: %.1f eV\nPer-pixel gain: %.1f eV'], Fr.fw, Fr.fwe, Fc.fw, Fc.fwe, Fr2.fw, Fc2.fw), ...
         'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 10);
    xlabel('Energy [eV]'); ylabel('Counts per 16 eV'); title('IMX290 ^{55}Fe Spectrum: Single Gain vs Per-Pixel Gain (10,000 pixels)');
    imx662_style(f);

    % 2. Mn K-alpha line (linear), 4 eV bins
    f = figure('Name', 'IMX290 Mn Ka line before and after gain correction', 'Position', [40 40 1000 650]);
    m = E.ec > 5550 & E.ec < 6250; mr = E.erw > 5550 & E.erw < 6250; sc = G;     % raw counts per 4 eV
    stairs(E.erw(mr), E.raw(mr)*sc, 'Color', [0.6 0.6 0.6], 'LineWidth', 1.4); hold on;
    stairs(E.ec(m), E.cor(m), 'Color', c.col.n1, 'LineWidth', 1.4);
    xx = linspace(5550, 6250, 500);
    plot(xx, gm(Fr.p, xx)*sc, ':', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.6);
    plot(xx, gm(Fc.p, xx), '--', 'Color', c.col.fit, 'LineWidth', 1.6);
    xline(c.E.MnKa, ':', 'Color', [.5 .5 .5]);
    legend({sprintf('Single gain: FWHM %.1f eV', Fr.fw), sprintf('Per-pixel gain: FWHM %.1f eV', Fc.fw), ...
            'Gaussian Fit (single gain)', 'Gaussian Fit (per-pixel gain)'}, 'Location', 'northwest');
    xlabel('Energy [eV]'); ylabel('Counts per 4 eV'); title('IMX290 Mn K\alpha Line');
    imx662_style(f);

    % 3. gain map + histogram
    f = figure('Name', 'IMX290 per-pixel gain map', 'Position', [40 40 1500 620]);
    tiledlayout(1, 2, 'TileSpacing', 'compact');
    nexttile; scatter(XY(:,1), XY(:,2), 9, g, 'filled'); axis image; set(gca, 'YDir', 'normal');
    colormap(gca, imx662_viridis(256)); cb = colorbar; cb.Label.String = 'Gain [ADU/eV]';
    caxis(prctile(g, [1 99])); xlabel('X Pixel'); ylabel('Y Pixel'); title('IMX290 Per-Pixel Gain (ECC), 10,000 Pixels');
    nexttile; [hc, xe] = histcounts(g, 0.288:0.0002:0.300); xc = xe(1:end-1) + 0.0001;
    bar(xc, hc, 1, 'FaceColor', c.col.bar, 'EdgeColor', 'none'); hold on;
    pg = imx662_gfit(xc, hc, [max(hc) G std(g)]);
    xx = linspace(0.288, 0.300, 500); plot(xx, gpk(pg, xx), '-', 'Color', c.col.fit, 'LineWidth', 2);
    xline(pg(2), '--', 'Color', c.col.reb, 'LineWidth', 1.5);
    legend({'Pixel Counts', 'Gaussian Fit', 'Gaussian Center'}, 'Location', 'northeast');
    text(0.03, 0.95, sprintf('Center = %.5f ADU/eV\nSigma = %.5f ADU/eV\nGain variation = %.3f %%\n\\sigma_{gain} = %.1f eV at 5900 eV', ...
         pg(2), pg(3), 100*pg(3)/pg(2), pg(3)/pg(2)*c.E.MnKa), 'Units', 'normalized', 'VerticalAlignment', 'top');
    xlabel('Gain [ADU/eV]'); ylabel('Number of Pixels'); title('IMX290 - Gain Histogram (ECC)');
    imx662_style(f);

    % 4. gain profiles
    f = figure('Name', 'IMX290 gain profiles', 'Position', [40 40 1300 520]);
    tiledlayout(1, 2, 'TileSpacing', 'compact');
    for k = 1:2
        nexttile; e = 0:50:1000; [~, ~, bi] = histcounts(XY(:,k), e); ok = bi > 0;
        mu = accumarray(bi(ok), g(ok)', [numel(e)-1 1], @mean, NaN);
        nn = accumarray(bi(ok), 1, [numel(e)-1 1]);
        se = accumarray(bi(ok), g(ok)', [numel(e)-1 1], @(v) std(v)/sqrt(numel(v)), NaN);
        keep = nn >= 20;
        errorbar(e(keep)+25, mu(keep), se(keep), 'o-', 'Color', c.col.n1, 'MarkerFaceColor', c.col.n1); hold on;
        yline(G, '--', 'Median gain', 'Color', c.col.fit);
        if k == 1, xlabel('X Pixel'); else, xlabel('Y Pixel'); end
        ylabel('Mean Gain [ADU/eV]'); title(sprintf('IMX290 Gain vs %s (50-pixel bins)', char('X'+k-1)));
        ylim(G + [-1 1]*0.0012);
    end
    imx662_style(f);

    % 5. two methods
    f = figure('Name', 'IMX290 ECC gain vs Mn Ka peak gain', 'Position', [40 40 800 700]);
    scatter(E.gPk, g, 6, c.col.n1, 'filled', 'MarkerFaceAlpha', 0.35); hold on; axis equal;
    lim = prctile([g E.gPk], [0.5 99.5]); plot(lim, lim, '--', 'Color', c.col.fit, 'LineWidth', 1.5);
    xlim(lim); ylim(lim); r = corrcoef(E.gPk, g);
    text(0.04, 0.95, sprintf('Pearson r = %.3f\nMedian difference = %.2f %%', r(1,2), 100*median((g - E.gPk)./E.gPk)), ...
         'Units', 'normalized', 'VerticalAlignment', 'top');
    xlabel('Gain from Mn K\alpha Peak Position [ADU/eV]'); ylabel('Gain from ECC [ADU/eV]');
    title('IMX290 Per-Pixel Gain: ECC vs Peak Position'); legend({'Pixels', '1:1'}, 'Location', 'southeast');
    imx662_style(f);

    % 6. accuracy vs events
    f = figure('Name', 'IMX290 gain accuracy vs events per pixel', 'Position', [40 40 1400 560]);
    tiledlayout(1, 2, 'TileSpacing', 'compact');
    nexttile; N = acc(:,1);
    loglog(N, acc(:,2), 'o-', 'Color', c.col.n1, 'MarkerFaceColor', c.col.n1, 'LineWidth', 1.4); hold on;
    nn = logspace(1, log10(1500), 200); expc = sqrt(1./nn - 1/E.nAll); expc = expc / expc(1) * acc(1,2);
    plot(nn, expc, '--', 'Color', c.col.fit, 'LineWidth', 1.3);
    yline(100*pg(3)/pg(2), ':', sprintf('Pixel-to-pixel gain variation (%.2f%%)', 100*pg(3)/pg(2)), 'Color', c.col.n1234, 'LineWidth', 1.6, ...
          'LabelHorizontalAlignment', 'right');
    legend({'Measured', 'Counting statistics \propto \surd(1/N - 1/N_{all})'}, 'Location', 'southwest');
    xlabel('Events per Pixel'); ylabel('Gain Error, Robust RMS [%]'); title('IMX290 Calibration Error vs Events per Pixel');
    xlim([8 1600]);
    nexttile; errorbar(N, acc(:,4), acc(:,5), 'o-', 'Color', c.col.n1, 'MarkerFaceColor', c.col.n1, 'LineWidth', 1.4); hold on;
    set(gca, 'XScale', 'log');
    yline(Fr.fw, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.3); yline(Fc.fw, '--', 'Color', c.col.fit, 'LineWidth', 1.3);
    yline(fanoMn, ':', 'Color', c.col.n1234, 'LineWidth', 1.6);
    legend({'Corrected with N-event gains', 'Single gain', sprintf('All events (~%d per pixel)', round(E.nAll)), 'Fano Limit'}, 'Location', 'northeast');
    xlabel('Events per Pixel'); ylabel('Mn K\alpha FWHM [eV]'); title('IMX290 Energy Resolution vs Events per Pixel');
    xlim([8 1600]); ylim([112 146]);
    imx662_style(f);

    % 7. ECC fit quality
    f = figure('Name', 'IMX290 ECC correlation coefficient', 'Position', [40 40 900 600]);
    histogram(E.rho, 200, 'FaceColor', c.col.bar, 'EdgeColor', 'none'); set(gca, 'YScale', 'log');
    xlabel('Correlation Coefficient (Pixel Spectrum vs Reference)'); ylabel('Number of Pixels');
    title('IMX290 ECC Fit Quality');
    imx662_style(f);

    % 8. comparison with Schneider et al. (2024)
    f = figure('Name', 'IMX290 Mn Ka resolution comparison', 'Position', [40 40 1100 600]);
    vals = [ben.fw1; Fr.fw Fc.fw];
    b = bar(vals, 'grouped'); b(1).FaceColor = [0.6 0.6 0.6]; b(2).FaceColor = c.col.n1; hold on;
    yline(fanoMn, ':', 'Fano Limit', 'Color', c.col.n1234, 'LineWidth', 1.6, 'LabelHorizontalAlignment', 'left');
    set(gca, 'XTickLabel', {'Schneider et al. (2024)', 'Carissma (2026)'});
    for i = 1:2, for j = 1:2, text(i + (j-1.5)*0.29, vals(i,j)+1.2, sprintf('%.1f', vals(i,j)), 'HorizontalAlignment', 'center'); end, end
    legend({'Single gain', 'Per-pixel gain'}, 'Location', 'northeast'); ylim([100 145]);
    ylabel('Mn K\alpha FWHM [eV]'); title('IMX290 Mn K\alpha Energy Resolution');
    imx662_style(f);
    imx662_save_figs('IMX290', OUT);
end
