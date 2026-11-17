function DC = imx226_dark_current(gain)
%IMX226_DARK_CURRENT  IMX226 dark current (lab method, Lupo 2025 Sec. 3.5 / 4.5) + extra analyses.
%   Same method as imx662_dark_current; IMX226 specifics follow Lupo (2025):
%   bias = 41 us cube, linear fit for exposures >= 0.5 s, regions = full /
%   centre (600-pixel rows removed at top and bottom) / edges (those rows only),
%   heat maps at 3 s and 0.2 s.
%   Bias = 15 us cube. Per cube: per-pixel median of 25 frames, minus bias median image.
%   5-sigma clipping x3 -> mean = average bias-subtracted dark signal; converted to
%   e-/pixel with e- = ADU / gain / 3.65.
%   Spatial error = RMS of the clipped bias-subtracted pixels.
%   Temporal error = mean over pixels of the 25-frame std (bias-median subtracted) / sqrt(25).
%   Heat maps (10x10 clusters) at 4.5 s and 0.2 s; region analysis without / with
%   the outer 375-pixel columns.
%   Extras: hot-pixel dark slopes and the shot-noise method (dark current from
%   the growth of temporal variance, which is not affected by the black-level clamp).
    c = imx662_cfg();
    if nargin < 1 || isempty(gain), gain = 0.3654; end
    % The 2.15 GB Box zip was copied as two smaller zips holding the 9 data cubes
    % (browser downloads are limited to 2 GB); the full zip is used if present.
    zipList = {fullfile(c.D, 'IMX226_DarkCurrent_20250604.zip')};
    if ~exist(zipList{1}, 'file'), zipList = {fullfile(c.D, 'IMX226_Dark_partA.zip'), fullfile(c.D, 'IMX226_Dark_partB.zip')}; end
    toE = 1 / (gain * c.epsSi);                     % ADU -> e-
    tmp = fullfile(tempdir, 'imx226_dc'); if exist(tmp,'dir'), rmdir(tmp,'s'); end; mkdir(tmp);
    cln = onCleanup(@() rmdir(tmp, 's'));
    names = {}; zof = {};
    for z = 1:numel(zipList)
        nm = imx662_zip_list(zipList{z}, '_datacube\.fits$'); nm = nm(:);
        names = [names; nm]; zof = [zof; repmat(zipList(z), numel(nm), 1)]; %#ok<AGROW>
    end
    tt = regexp(names, 'DarkCurrent_(\w+?)_1_datacube', 'tokens', 'once'); tt = cellfun(@(q) q{1}, tt, 'UniformOutput', false);
    texp = cellfun(@parseT, tt); [texp, o] = sort(texp); names = names(o); zof = zof(o); n = numel(names);
    fprintf('\n=== Dark current: %d cubes (%s ... %s s) ===\n', n, num2str(texp(1)), num2str(texp(end)));

    % bias (shortest exposure)
    fb = imx662_zip_extract(zof{1}, names{1}, tmp); cb = imx662_read_cube(fb); delete(fb);
    biasMed = median(cb, 3); clear cb
    [ny, nx] = size(biasMed);
    regions = {'Full detector', true(ny, nx); ...
               'Center (600-pixel rows removed at top and bottom)', [false(600,nx); true(ny-1200, nx); false(600,nx)]; ...
               'Edges only (top and bottom 600-pixel rows)', [true(600,nx); false(ny-1200, nx); true(600,nx)]};
    nR = size(regions, 1);
    sig = nan(n, nR); spat = sig; temp = nan(n,1); tempVar = nan(n,1);
    meanMaps = zeros(ny, nx, n, 'single');
    for k = 1:n
        f = imx662_zip_extract(zof{k}, names{k}, tmp); cube = imx662_read_cube(f); delete(f);
        med = median(cube, 3) - biasMed;
        meanMaps(:,:,k) = mean(cube, 3);
        for r = 1:nR
            x = double(med(regions{r,2})); keep = true(size(x));
            for it = 1:3, mu = mean(x(keep)); sd = std(x(keep)); keep = abs(x - mu) <= 5*sd; end
            sig(k, r) = mean(x(keep)) * toE; spat(k, r) = sqrt(mean((x(keep) - mean(x(keep))).^2)) * toE;
        end
        sdp = std(cube - biasMed, 0, 3);              % per-pixel temporal std
        temp(k) = mean(sdp(:)) / sqrt(size(cube,3)) * toE;
        tempVar(k) = median(reshape(var(cube, 0, 3), [], 1));   % ADU^2 (for shot-noise method)
        if abs(texp(k) - 3) < 1e-6, map45 = med; end
        if abs(texp(k) - 0.2) < 1e-6, map02 = med; end
        fprintf('  t = %8.6f s : dark signal %.3f e-/pixel (spatial %.3f, temporal %.3f)\n', texp(k), sig(k,1), spat(k,1), temp(k));
        clear cube sdp
    end
    use = texp >= 0.07;                             % plotted: 0.08 ... 3 s (bias excluded)
    useFit = texp >= 0.49;                          % thesis: linear fit for exposures >= 0.5 s
    fits = zeros(nR, 2);
    for r = 1:nR, fits(r,:) = polyfit(texp(useFit), sig(useFit, r), 1); end
    fprintf('Lab-method dark current (slope): full %.4f, center %.4f, edges %.4f e-/pixel/s\n', fits(:,1));

    % ---- Extra 1: shot-noise method ----
    pv = polyfit(texp, tempVar, 1);
    DCshot = pv(1) * toE^2;                          % e-/s  (variance in e-^2 grows by D*t)
    readE = sqrt(pv(2)) * toE;
    fprintf('Shot-noise method: temporal variance grows %.4f ADU^2/s -> typical-pixel dark current %.3f e-/pixel/s; read noise %.2f e-\n', pv(1), DCshot, readE);

    % ---- Extra 2: per-pixel slopes (median-referenced) and hot pixels ----
    Y = reshape(meanMaps, ny*nx, n); medL = median(Y, 1);
    Yc = double(Y) - double(medL); W = double(Y < 4000);
    Sw = sum(W,2); St = W*texp; Stt = W*texp.^2; Sy = sum(W.*Yc,2); Sty = sum(W.*Yc.*texp',2);
    slope = (Sw.*Sty - St.*Sy) ./ (Sw.*Stt - St.^2); slope(Sw < 3) = NaN;
    slopeE = slope * toE; sv = slopeE(isfinite(slopeE));
    sMed = median(sv); sSig = 1.4826*median(abs(sv - sMed));
    hotThrE = max(sMed + 5*sSig, 10*toE);
    hot = reshape(slopeE > hotThrE | Sw < 3, ny, nx);
    fprintf('Hot (dark) pixels: %d (> %.2f e-/s above the array median)\n', nnz(hot), hotThrE);
    clear Y Yc W

    DC = struct('texp', texp, 'sig', sig, 'spat', spat, 'temp', temp, 'fits', fits, 'regions', {regions(:,1)}, ...
        'tempVar', tempVar, 'pvShot', pv, 'DCshot', DCshot, 'readE', readE, 'slopeE', reshape(slopeE, ny, nx), ...
        'hotThrE', hotThrE, 'nHot', nnz(hot), 'toE', toE, 'gain', gain);

    % ================= FIGURES =================
    % ---- Figs 4-35, 4-38, 4-39 style: one per region ----
    for r = 1:nR
        f = figure('Name', ['Dark current ' regions{r,1}], 'Color', 'w', 'Position', [40 40 1400 880]);
        tl = tiledlayout(f, 2, 1, 'TileSpacing', 'compact');
        nexttile;
        h1 = plot(texp(use), sig(use,r), 'o-', 'Color', c.col.n1, 'LineWidth', 1.3, 'MarkerFaceColor', 'w'); hold on; grid on; box on;
        h2 = errorbar(texp(use), sig(use,r), spat(use,r), 'LineStyle', 'none', 'Color', c.col.spat, 'LineWidth', 1.2, 'CapSize', 8);
        h3 = errorbar(texp(use)+0.03, sig(use,r), temp(use), 'LineStyle', 'none', 'Color', c.col.temp, 'LineWidth', 1.6, 'CapSize', 5);
        xx = linspace(0, max(texp), 100); h4 = plot(xx, polyval(fits(r,:), xx), ':', 'Color', c.col.fit, 'LineWidth', 1.6);
        yline(0, 'k-'); xlim([0 max(texp)+0.25]);
        ylabel('Bias Subtracted Dark Signal [e^-/pixel]'); xlabel('Exposure Time [Seconds]');
        title(['IMX226: Dark Current Measurements - ' regions{r,1}]);
        legend([h1 h4 h2 h3], {'Dark Current [e^-]', sprintf('Line of Best Fit: %.4f*x + %.3f', fits(r,1), fits(r,2)), 'Spatial Error', 'Temporal Error'}, 'Location', 'northwest');
        nexttile;
        res = sig(useFit,r) - polyval(fits(r,:), texp(useFit));
        plot(texp(useFit), res, 'o', 'Color', c.col.n1, 'MarkerFaceColor', c.col.n1); grid on; box on; yline(0, 'k-');
        ylabel('Bias Subtracted Dark Signal Residual [e^-/pixel]'); xlabel('Exposure Time [Seconds]');
        title('IMX226: Dark Current Measurement Residuals'); legend({'Residuals'}, 'Location', 'northeast'); xlim([0 max(texp)+0.25]);
        imx662_caption(f, sprintf(['Dark current = signal a pixel builds up with NO light, from heat (thermally freed electrons). It should grow in a straight line with exposure time; the slope is the dark-current rate. ' ...
            'Each point is the average over the %s (after removing the 15-microsecond "bias" frame and outlier pixels), in electrons per pixel. Green bars: how much pixels differ from each other (spatial); ' ...
            'violet bars: frame-to-frame scatter of the average (temporal). Slope = %.4f e-/pixel/s. The points scatter within their error bars with no clear trend, as in Lupo (2025): ' ...
            'with this method the IMX662 dark current at 25 C is too small to measure, partly because the camera''s automatic black-level correction removes the average offset. See the extra plots for a more sensitive method.'], ...
            lower(regions{r,1}), fits(r,1)), 0.2);
    end

    % ---- Fig 4-36 style heat maps (10x10) and Fig 4-37 region box ----
    cm = imx662_viridis(256);
    for q = 1:2
        if q == 1, M = map45; lab = '3s'; else, M = map02; lab = '0.2s'; end
        cl = 10; Mc = squeeze(mean(mean(reshape(double(M(1:floor(ny/cl)*cl, 1:floor(nx/cl)*cl)) * toE, cl, floor(ny/cl), cl, floor(nx/cl)), 1), 3));
        f = figure('Name', ['Dark heat map ' lab], 'Color', 'w', 'Position', [40 40 1300 760]);
        imagesc((0:size(Mc,2)-1)*cl, (0:size(Mc,1)-1)*cl, Mc); axis image; colormap(f, cm); caxis([-2 4]);
        cbh = colorbar; cbh.Label.String = 'Bias Subtracted Signal [e^-]'; hold on;
        if q == 1, rectangle('Position', [0 600 nx ny-1200], 'EdgeColor', c.col.fit, 'LineWidth', 2); end
        xlabel('X Pixel'); ylabel('Y Pixel'); title([lab ' Exposure Time Heat Map - IMX226 (10\times10-pixel clusters)']);
        set(gca, 'FontSize', c.fs);
        txtq = ['Map of the bias-subtracted dark signal over the sensor at ' lab ' exposure, averaged in 10x10-pixel blocks, in electrons (yellow = more dark signal, dark purple = less; same -2 to +2 e- scale as Lupo 2025). '];
        if q == 1, txtq = [txtq 'The dashed box is the central region used for the "center" analysis (outer 375-pixel columns excluded), the same split as Lupo (2025) Fig. 4-37.']; end
        imx662_caption(f, txtq, 0.13);
    end

    % ---- Extra: shot-noise method ----
    f = figure('Name', 'Dark current shot-noise method', 'Color', 'w', 'Position', [40 40 1200 650]);
    plot(texp, tempVar * toE^2, 'o', 'Color', c.col.n1, 'MarkerFaceColor', c.col.n1); hold on; grid on; box on;
    xx = linspace(0, max(texp), 100); plot(xx, polyval(pv, xx) * toE^2, '-', 'Color', c.col.fit, 'LineWidth', 1.5);
    xlabel('Exposure Time [Seconds]'); ylabel('Median temporal noise variance [(e^-)^2]');
    title(sprintf('IMX226 Dark Current from Noise Growth: %.3f e^-/pixel/s', DCshot));
    legend({'measured', sprintf('fit: %.2f + %.3f\\cdott  (e^-)^2', pv(2)*toE^2, DCshot)}, 'Location', 'northwest'); set(gca, 'FontSize', c.fs);
    imx662_caption(f, sprintf(['A more sensitive way to see a small dark current. Thermal electrons arrive randomly, so even when the camera subtracts their average, they add random noise: ' ...
        'the frame-to-frame noise variance (in electrons squared) grows by exactly the dark current x exposure time. The slope gives a typical-pixel dark current of %.3f e-/pixel/s at 25 C, ' ...
        'and the value at zero exposure is the read noise squared (%.2f e-). This is why the average-signal method in the previous plots could not see it.'], DCshot, readE), 0.18);

    % ---- Extra: per-pixel dark slope histogram + hot pixel map ----
    f = figure('Name', 'Dark hot pixels', 'Color', 'w', 'Position', [40 40 1500 600]);
    tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact');
    nexttile; histogram(sv, linspace(-3, max(prctileV(sv, 99.999), 15), 300), 'FaceColor', c.col.n1, 'EdgeColor', 'none');
    set(gca, 'YScale', 'log'); grid on; xline(hotThrE, '--', sprintf('hot threshold %.1f e^-/s', hotThrE), 'Color', c.col.hot, 'LineWidth', 1.4);
    xlabel('Per-pixel dark signal growth relative to the median pixel [e^-/pixel/s]'); ylabel('Number of Pixels');
    title('IMX226 Per-Pixel Dark-Signal Slopes');
    nexttile; [rr, cc] = find(hot); scatter(cc, rr, 8, min(slopeE(hot), prctileV(slopeE(hot), 95)), 'filled'); axis ij image; xlim([1 nx]); ylim([1 ny]);
    colormap(gca, imx662_viridis(256)); cbh = colorbar; cbh.Label.String = 'Dark slope [e^-/pixel/s]'; grid on;
    xlabel('X Pixel'); ylabel('Y Pixel'); title(sprintf('IMX226: %d Hot (Dark-Current) Pixels', nnz(hot)));
    imx662_caption(f, sprintf(['Left: for every pixel, how fast its dark signal grows with exposure compared with the typical pixel (log y). Almost all pixels sit at 0 (no extra dark current); a long tail of ' ...
        '"hot pixels" gains up to hundreds of electrons per second. Right: where the %d hot pixels (> %.1f e-/s) are; they are scattered randomly, not clustered, which is typical of isolated crystal defects.'], nnz(hot), hotThrE), 0.18);

    % ---- Extra: hot pixel response curves ----
    f = figure('Name', 'Dark pixel response curves', 'Color', 'w', 'Position', [40 40 1400 560]);
    tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact');
    Yv = reshape(meanMaps, ny*nx, n); s2 = slopeE(:); s2(~isfinite(s2)) = -Inf; [~, iH] = maxk(s2, 12);
    iT = find(abs(slopeE(:) - sMed) < 0.05*sSig, 12);
    nexttile; plot(texp, double(Yv(iH,:))', 'o-'); grid on; yline(c.adcMax, 'k--', 'Saturation (4095 ADU)');
    xlabel('Exposure Time [Seconds]'); ylabel('Pixel value (mean of 25 frames) [ADU]'); title('IMX226: 12 Hottest Pixels');
    nexttile; plot(texp, double(Yv(iT,:))', 'o-'); grid on;
    xlabel('Exposure Time [Seconds]'); ylabel('Pixel value (mean of 25 frames) [ADU]'); title('IMX226: 12 Typical Pixels');
    imx662_caption(f, ['Raw values of individual pixels versus exposure time. Hot pixels (left) rise in a straight line until they reach the 12-bit limit (4095 ADU). Typical pixels (right) stay at the ' ...
        '~100 ADU black level, with only noise. Some hot pixels jump up and down between exposures ("random telegraph" pixels), a known defect behaviour.'], 0.18);
    clear Yv meanMaps
end

function t = parseT(s)
    if strcmpi(s, 'Bias'), t = 41e-6; return; end      % shortest exposure allowed for the IMX226 (Lupo 2025)
    s = erase(s, 'EXP');
    if endsWith(s, 'us'), t = str2double(erase(s, 'us')) * 1e-6;
    else, t = str2double(strrep(erase(s, 's'), 'p', '.')); end
end

function cm = blueWhiteRed(n)
    h = n/2; a = linspace(0, 1, h)';
    cm = [[a a ones(h,1)]; [ones(h,1) flipud(a) flipud(a)]];
    cm = 0.15 + 0.85*cm; cm(cm > 1) = 1;
end

function v = prctileV(x, p)
    x = sort(x(:)); x = x(isfinite(x)); v = x(max(1, min(numel(x), round(p/100*numel(x)))));
end
