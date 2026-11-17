function imx662_plot_xray(R)
%IMX662_PLOT_XRAY  All X-ray response figures (lab standard + extras) from imx662_xray_events results.
    c = imx662_cfg(); isFe = R.isFe;
    % sensor-specific settings (defaults = IMX662; IMX226 results carry their own)
    sn = 'IMX662'; if isfield(R, 'sensor'), sn = R.sensor; end
    if isfield(R, 'cfg'), fn = fieldnames(R.cfg); for i = 1:numel(fn), c.(fn{i}) = R.cfg.(fn{i}); end, end
    if ~isfield(c, 'cluster'), c.cluster = 3; end
    xMaxEV = R.xMaxEV; xc = R.xc; hc = R.hc;
    pn = @(s) strrep(strrep(s, '\alpha', '-alpha'), '\beta', '-beta');   % plain names for captions
    L = R.L; S = R.S; g = R.gain; ev = (1:numel(S.n1))' / g; src = R.src;
    if isFe, srcTxt = 'Fe-55'; else, srcTxt = 'Po-210'; end
    gm = @(p, x) p(1)*exp(-0.5*((x-p(2))/p(3)).^2);

    % ---- Pixel event-count histogram / hot-pixel threshold ----
    f = figure('Name', [srcTxt ' pixel count histogram'], 'Color', 'w', 'Position', [40 40 950 650]);
    b = bar(xc, hc, 1, 'FaceColor', c.col.bar, 'EdgeColor', 'none', 'FaceAlpha', 0.75); hold on; grid on; box on;
    if isFe
        xx = linspace(0, max(R.hotThr*1.6, R.pg(2)*2.5), 800); plot(xx, gm(R.pg, xx), '-', 'Color', c.col.fit, 'LineWidth', 2);
        xline(R.pg(2), '--', 'Color', c.col.reb, 'LineWidth', 1.5);
        xline(R.hotThr, ':', 'Color', c.col.hot, 'LineWidth', 2);
        legend({'Counts', 'Gaussian Fit', 'Gaussian Center', sprintf('Hot-pixel threshold = center + 3\\sigma = %.0f', R.hotThr)}, 'Location', 'northeast');
        text(0.02, 0.95, sprintf('Center = %.2f \\pm %.2f\nSigma = %.2f', R.pg(2), R.pge(2), R.pg(3)), 'Units', 'normalized', ...
             'VerticalAlignment', 'top', 'BackgroundColor', 'w', 'EdgeColor', [.6 .6 .6]);
        xlim([0 max(R.hotThr*1.6, R.pg(2)*2.5)]);
        cap = sprintf(['Every pixel on the sensor was counted for how many X-ray events it recorded over the whole %.1f-hour Fe-55 run. ' ...
               'Most pixels follow a bell curve (Gaussian) because the source lights the sensor fairly evenly. Pixels far to the right of the curve ' ...
               '(more than 3 standard deviations above its center, i.e. > %.0f events) fire too often to be real X-rays: these are "hot pixels" and are removed ' ...
               'from all spectra. %d pixels (%.2f%% of the sensor) were flagged.'], R.tTot/3600, R.hotThr, R.nHot, R.hotPct);
    else
        set(gca, 'YScale', 'log'); xline(R.hotThr, ':', 'Color', c.col.hot, 'LineWidth', 2);
        legend({'Counts', 'Hot-pixel threshold (\geq 10 events)'}, 'Location', 'northeast'); xlim([-0.5 max(40, R.hotThr*3)]);
        cap = sprintf(['Number of events recorded by each pixel over the %.1f-hour Po-210 run (log scale). The Po-210 source is weak, so a real pixel sees only ' ...
               'a few X-rays; any pixel with 10 or more events is treated as a hot pixel (lab rule) and removed. %d pixels (%.2f%%) were flagged.'], ...
               R.tTot/3600, R.nHot, R.hotPct);
    end
    xlabel('Counts per Pixel'); ylabel('Number of Pixels'); title([sn ' ' srcTxt ' Pixel Count']);
    set(gca, 'FontSize', c.fs); imx662_caption(f, cap, 0.2);

    % ---- Spectrum (lab-standard plot) ----
    f = figure('Name', [srcTxt ' spectrum'], 'Color', 'w', 'Position', [40 40 1500 820]);
    sm = @(h) max(h, 0.5);
    hG = stairs(ev, sm(S.n1234), 'Color', c.col.n1234, 'LineWidth', 2.2); hold on;
    hB = stairs(ev, sm(S.n1), 'Color', c.col.n1, 'LineWidth', 1.1);
    hA = stairs(ev, sm(S.all), 'Color', c.col.all, 'LineWidth', 0.8);
    set(gca, 'YScale', 'log'); grid on; box on;
    yl = [1 max(S.all)*10^2.4]; ylim(yl); xlim([0 xMaxEV]);
    for j = 1:numel(L)
        for cls = {'n1','n1234'}
            q = L(j).(cls{1}); if ~isfinite(q.mu) || q.A <= 0, continue; end
            xx = linspace(q.mu - 3*q.sig, q.mu + 3*q.sig, 200);
            if strcmp(cls{1}, 'n1'), cc = c.col.n1; else, cc = c.col.n1234; end
            plot(xx / g, max(gm([q.A q.mu q.sig], xx) + max(q.bkg, 0), 0.5), '--', 'Color', cc*0.75, 'LineWidth', 1.4);
        end
        xline(L(j).n1.mu / g, ':', 'Color', [.5 .5 .5]);
        q1 = L(j).n1; q4 = L(j).n1234;
        txt = sprintf('%s line\nn1:\nCenter = %.1f eV\nSigma = %.1f eV\nFWHM = %.1f \\pm %.1f eV\nCounts = %.0f\n\nn1234:\nCenter = %.1f eV\nSigma = %.1f eV\nFWHM = %.1f \\pm %.1f eV\nCounts = %.0f', ...
            L(j).short, q1.centerEV, q1.sigmaEV, q1.fwhmEV, q1.fwhmEVerr, q1.area, q4.centerEV, q4.sigmaEV, q4.fwhmEV, q4.fwhmEVerr, q4.area);
        xt = 0.02 + (j-1) * 0.98/numel(L);
        text(xt, 0.985, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8.5, 'Color', [0.1 0.1 0.1]);
    end
    if isFe
        xline([c.E.Ti c.E.MnKbEsc], '-.', {'Ti K\alpha?', 'Si esc (Mn K\beta)'}, 'Color', [.6 .6 .6], 'LabelVerticalAlignment', 'bottom');
    end
    legend([hA hB hG], {'all', 'n1', 'n1234'}, 'Location', 'southwest');
    xlabel('Energy [eV]'); ylabel('Counts'); title([sn ' ' srcTxt ' Spectrum']);
    set(gca, 'FontSize', c.fs);
    if isFe
        cap = sprintf(['Energy spectrum of all X-rays detected from the Fe-55 source (hot pixels removed). Each detected X-ray is placed by its energy (x, in eV; converted from ADU with the gain %.4f ADU/eV ' ...
          'measured from the Mn K-alpha peak) and the y-axis counts how many X-rays had that energy (log scale). Peaks are emission lines: Mn K-alpha (5900 eV) and Mn K-beta (6490 eV) from the source, ' ...
          'Si K-alpha (1740 eV) from fluorescence of the silicon sensor, and the Si "escape" peak (~4160 eV) where a Si X-ray left the sensor and carried away 1740 eV. ' ...
          'Blue (n1) = X-rays whose charge landed in a single pixel, green (n1234) = charge in 1-4 pixels, black (all) = up to 9 pixels. Dashed curves are Gaussian fits; FWHM (peak width) is the energy resolution, smaller is better.'], g);
    else
        cap = sprintf(['Energy spectrum from the Po-210 source (hot pixels removed). Po-210 alpha particles hit a Teflon (C, F) target, which emits low-energy X-rays: Carbon K-alpha (277 eV) and Fluorine K-alpha (677 eV); ' ...
          'a weak Oxygen K-alpha bump (~526 eV) comes from contamination. Energies use the Fe-55 gain (%.4f ADU/eV). The sharp cut near %.0f eV is the 30 ADU event-detection threshold. ' ...
          'Blue = single-pixel (n1), green = 1-4 pixels (n1234), black = up to 9 pixels (all). Dashed curves are Gaussian fits; FWHM = energy resolution.'], g, 30/g);
    end
    imx662_caption(f, cap, 0.17);

    % ---- Zoomed linear views of each fitted line (extra, for clarity) ----
    f = figure('Name', [srcTxt ' line fits zoom'], 'Color', 'w', 'Position', [40 40 1500 520]);
    tl = tiledlayout(f, 1, numel(L), 'TileSpacing', 'compact');
    for j = 1:numel(L)
        nexttile; q = L(j).n1; xr = [q.mu - 5*q.sig, q.mu + 5*q.sig];
        sel = (1:numel(S.n1))' >= xr(1) & (1:numel(S.n1))' <= xr(2);
        stairs(ev(sel), S.n1(sel), 'Color', c.col.n1, 'LineWidth', 1.2); hold on; grid on;
        stairs(ev(sel), S.n1234(sel), 'Color', c.col.n1234, 'LineWidth', 1.2);
        xx = linspace(xr(1), xr(2), 300);
        plot(xx/g, gm([q.A q.mu q.sig], xx) + q.bkg, '-', 'Color', c.col.fit, 'LineWidth', 1.5);
        xline(L(j).E, '--k');
        ax = gca; ax.YAxis.Exponent = 0; ax.XAxis.Exponent = 0;
        title(sprintf('%s', L(j).name)); subtitle(sprintf('FWHM %.1f eV (n1)', q.fwhmEV)); xlabel('Energy [eV]'); ylabel('Counts');
    end
    lg = legend({'n1', 'n1234', 'Gaussian fit (n1)', 'Expected line energy'}, 'Orientation', 'horizontal'); lg.Layout.Tile = 'south';
    imx662_caption(f, ['Close-up of each emission line on a linear scale. The red curve is the Gaussian fit to the single-pixel peak; the dashed black line marks the ' ...
        'textbook energy of that line. If the peak sits on the dashed line, the energy calibration is correct there. The width of the peak (FWHM) is the energy resolution.'], 0.2);

    % ---- Hot-pixel spectrum ----
    f = figure('Name', [srcTxt ' hot pixel spectrum'], 'Color', 'w', 'Position', [40 40 1300 650]);
    Sh = S.hot; Sh(end) = 0;
    stairs(ev, max(Sh, 0.5), 'Color', c.col.hot, 'LineWidth', 1.2); set(gca, 'YScale', 'log'); grid on; box on; ylim([0.5 max(Sh)*3]);
    xlim([0 min(xMaxEV * (1 + 0.8*isFe), ev(end-1))]); xlabel('Energy [eV]'); ylabel('Counts');
    title([sn ' ' srcTxt ' Spectrum of Hot Pixels']);
    for j = 1:numel(L), xline(L(j).E, ':', L(j).name, 'Color', [.4 .4 .4]); end
    if isFe, xline(2*[c.E.MnKa c.E.MnKb], '-.', {'2 \times Mn K\alpha (pile-up)', '2 \times Mn K\beta'}, 'Color', [.5 .5 .5]); end
    set(gca, 'FontSize', c.fs);
    imx662_caption(f, ['Spectrum of only the events from the removed hot pixels. Most counts pile up at the lowest energies (electronic noise mistaken for signal), which is why these ' ...
        'pixels are removed. Some real emission-line X-rays also appear because the hot-pixel cut is not perfect; for Fe-55, bumps near twice the Mn energies are "pile-up", ' ...
        'where two X-rays hit the same pixel in one frame and are recorded as one event with double energy.'], 0.18);

    % ---- Light curves ----
    for j = find([L.lc])
        rate = L(j).lcAll(R.present) / c.frameTime / c.areaMM2;
        t = (0:numel(rate)-1)' * c.frameTime;
        nb = max(round(numel(rate)*0.001), 1);  nbin = floor(numel(rate)/nb);
        rb = mean(reshape(rate(1:nbin*nb), nb, nbin), 1)'; tb = mean(reshape(t(1:nbin*nb), nb, nbin), 1)';
        f = figure('Name', [srcTxt ' light curve ' L(j).short], 'Color', 'w', 'Position', [40 40 1500 560]);
        tl = tiledlayout(f, 1, 4, 'TileSpacing', 'compact');
        ax1 = nexttile([1 3]);
        plot(t, rate, '-', 'Color', c.col.raw, 'LineWidth', 0.5); hold on; grid on; box on;
        plot(tb, rb, '-', 'Color', c.col.reb, 'LineWidth', 1.6);
        yline(mean(rate), '--', sprintf('mean %.4g', mean(rate)), 'Color', c.col.reb, 'LabelHorizontalAlignment', 'left', 'LineWidth', 1.2);
        xlabel('Seconds'); ylabel('Count Rate for the Detector [Counts/Sec/mm^2]'); ax1.XAxis.Exponent = 0;
        title([sn ' Light Curve for ' L(j).name]);
        legend({[L(j).short ' Curve'], ['Re-Binned ' L(j).short ' Curve']}, 'Location', 'southoutside', 'Orientation', 'horizontal');
        ax2 = nexttile;
        cts = L(j).lcAll(R.present); k2r = 1 / c.frameTime / c.areaMM2;
        if max(cts) <= 60
            histogram(rate, ((0:max(cts)+1) - 0.5) * k2r, 'Orientation', 'horizontal', 'FaceColor', c.col.hist, 'EdgeColor', 'w');
        else
            histogram(rate, 40, 'Orientation', 'horizontal', 'FaceColor', c.col.hist, 'EdgeColor', 'none');
        end
        grid on; ax2 = gca; ax2.XAxis.Exponent = 0;
        xlabel('Number of Occurrences'); linkaxes([ax1 ax2], 'y'); ax2.YTickLabel = [];
        set([ax1 ax2], 'FontSize', c.fs);
        imx662_caption(f, sprintf(['"Light curve" = how many %s X-rays (energies within 3 sigma of the line center, all event sizes, hot pixels removed) the sensor detected per second per mm^2, ' ...
            'versus time since the start of the run. Red: every frame; blue: averaged into bins of 0.1%% of the run to show the trend. A flat curve means the source and detector were stable. ' ...
            'Right: histogram of the per-frame rates; its width is mostly counting (Poisson) noise. Mean rate %.4g counts/s/mm^2. ' ...
            'For weak lines only 0, 1, 2 ... X-rays arrive per frame, so the raw rate takes discrete steps of %.3g counts/s/mm^2 (one X-ray per frame); the blue curve is the meaningful average.'], ...
            pn(L(j).name), mean(rate), 1/c.frameTime/c.areaMM2), 0.2);
    end

    % ---- Count-rate heat maps (3x3 clusters) ----
    for j = find([L.lc])
        cl = c.cluster; ny3 = floor(c.ny/cl); nx3 = floor(c.nx/cl);
        if isfield(L(j), 'mapC') && ~isempty(L(j).mapC)          % already summed in cl x cl clusters
            sC = L(j).mapC; nG = cl^2 - R.hotC;
        else
            m = reshape(L(j).map, c.ny, c.nx); h = R.hotMap;
            mm = m(1:ny3*cl, 1:nx3*cl); hh = ~h(1:ny3*cl, 1:nx3*cl);
            sC = squeeze(sum(sum(reshape(mm, cl, ny3, cl, nx3), 1), 3));
            nG = squeeze(sum(sum(reshape(hh, cl, ny3, cl, nx3), 1), 3));
        end
        rateMap = sC ./ (R.tTot * max(nG,1) * c.pixMM^2); rateMap(nG == 0) = NaN;
        f = figure('Name', [srcTxt ' heat map ' L(j).short], 'Color', 'w', 'Position', [40 40 1300 760]);
        imagesc((0:nx3-1)*cl, (0:ny3-1)*cl, rateMap); axis xy image; colormap(f, imx662_viridis(256)); cb = colorbar;
        cb.Label.String = 'Count Rate [Counts/Second/mm^2]';
        lo = prctile0(rateMap(isfinite(rateMap)), 0.5); hi = prctile0(rateMap(isfinite(rateMap)), 99.5);
        if hi > lo, caxis([lo hi]); end
        xlabel('X pixel'); ylabel('Y pixel'); title([sn ' Count Rate Map for ' L(j).name sprintf(' (%d\\times%d-pixel clusters)', cl, cl)]);
        set(gca, 'FontSize', c.fs);
        imx662_caption(f, sprintf(['Where on the sensor the %s X-rays landed. The sensor was divided into 3x3-pixel blocks (75.7 square micrometers each) and the colour is the count rate ' ...
            'in that block, per second per mm^2. Smooth gradients come from the source geometry (the source is closer to / tilted toward one side); isolated dark spots would indicate dead regions or ' ...
            'surface contamination; white = blocks made only of hot pixels.'], pn(L(j).name)), 0.15);
    end

    % ---- Event size distribution (extra) ----
    f = figure('Name', [srcTxt ' event size'], 'Color', 'w', 'Position', [40 40 900 600]);
    bar(1:9, 100*R.npixHist/sum(R.npixHist), 'FaceColor', c.col.n1); grid on;
    for k = 1:9, text(k, 100*R.npixHist(k)/sum(R.npixHist), sprintf('%.1f%%', 100*R.npixHist(k)/sum(R.npixHist)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom'); end
    xlabel('Number of pixels in the event (above the 9 ADU split threshold)'); ylabel('Percent of events [%]');
    title([sn ' ' srcTxt ': Event Size Distribution']); set(gca, 'FontSize', c.fs);
    imx662_caption(f, ['When an X-ray is absorbed it makes a small cloud of electrons. If the cloud stays in one pixel the event is "single-pixel" (best energy resolution); ' ...
        'if it spreads into neighbours it is a "multi-pixel" event. A high single-pixel fraction means the pixel walls confine charge well (STARVIS 2 design).'], 0.2);

    % ---- Gain variance (Fe-55) ----
    if isFe && isfield(R.GV, 'center')
        GV = R.GV; bc = (GV.be(1:end-1) + GV.be(2:end))/2;
        f = figure('Name', 'Fe-55 gain variance', 'Color', 'w', 'Position', [40 40 950 650]);
        bar(bc, GV.hg, 1, 'FaceColor', c.col.bar, 'EdgeColor', 'none', 'FaceAlpha', 0.75); hold on; grid on; box on;
        xx = linspace(GV.be(1), GV.be(end), 400); plot(xx, gm(GV.fit, xx), '-', 'Color', c.col.fit, 'LineWidth', 2);
        xline(GV.center, '--', 'Color', c.col.reb, 'LineWidth', 1.5); xline(R.gain, ':', 'Color', [0 0 0], 'LineWidth', 1.5);
        legend({'Pixel Counts', 'Gaussian Fit', 'Gaussian Center', sprintf('Detector gain %.4f ADU/eV', R.gain)}, 'Location', 'northeast');
        text(0.02, 0.95, sprintf('Center = %.5f\nSigma = %.5f\nGain variation = %.2f%%\n\\sigma_{gain} = %.1f eV at 5900 eV\nPixels used: %d', ...
             GV.center, GV.sigma, GV.pct, GV.sigmaEV, GV.nPix), 'Units', 'normalized', 'VerticalAlignment', 'top', 'BackgroundColor', 'w', 'EdgeColor', [.6 .6 .6]);
        xlabel('Gain [ADU/eV]'); ylabel('Number of Pixels'); title([sn ' - Gain Histogram of Mn K\alpha']); set(gca, 'FontSize', c.fs);
        imx662_caption(f, sprintf(['Each pixel converts collected charge to a digital number (ADU) with slightly different "gain". For every pixel with more than 100 single-pixel Mn K-alpha X-rays ' ...
            '(%d pixels), its gain = its average Mn K-alpha signal [ADU] / 5900 eV. The spread of this histogram (sigma/center = %.2f%%) is the pixel-to-pixel gain variation; ' ...
            'multiplied by 5900 eV it gives %.1f eV of extra line broadening when one gain is used for the whole sensor.'], GV.nPix, GV.pct, GV.sigmaEV), 0.2);
    end
end

function v = prctile0(x, p)
    x = sort(x(:)); v = x(max(1, min(numel(x), round(p/100*numel(x)))));
end
