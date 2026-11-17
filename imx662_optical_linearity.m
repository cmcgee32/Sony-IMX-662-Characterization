function OL = imx662_optical_linearity(gain)
%IMX662_OPTICAL_LINEARITY  Optical linearity (lab method, Lupo 2025 Sec. 3.4.2) + photon transfer.
%   100 pairs of 25-frame cubes, shutter closed (dark) and open (integrating
%   sphere), exposure from 2.0 s down to 0.02 s in 0.02 s steps (+16 us pair).
%   Per cube: per-pixel median over the 25 frames. Signal = open - closed.
%   5-sigma clipping repeated 3 times, then mean over pixels = mean dark-subtracted signal.
%   Fits: 5-95% of saturation (Fig 4-19), < 2000 ADU (Fig 4-20), 2000-4000 ADU (Fig 4-21).
%   Linearity error (Eq. 3.4): E = (max(r)-min(r)) / (max(l)-min(l)) * 100.
%   Extra: photon-transfer curve from frame-to-frame variance -> conversion gain in e-/ADU.
    c = imx662_cfg();
    if nargin < 1 || isempty(gain), gain = c.gainRef; end
    tmp = fullfile(tempdir, 'imx662_ol'); if exist(tmp,'dir'), rmdir(tmp,'s'); end; mkdir(tmp);
    cln = onCleanup(@() rmdir(tmp, 's'));
    names = imx662_zip_list(c.zipLin, 'datacube_Shutter(Open|Closed)_[\d.]+\.fits$');
    tok = regexp(names, 'datacube_Shutter(Open|Closed)_([\d.]+)\.fits$', 'tokens', 'once');
    isOpen = cellfun(@(q) strcmp(q{1}, 'Open'), tok); tus = cellfun(@(q) str2double(q{2}), tok);
    tex = unique(tus); n = numel(tex);
    fprintf('\n=== Optical linearity: %d exposure pairs ===\n', n);
    S = nan(n,1); Sstd = S; Vshot = S; Smed = S; dark = S; satPix = S; t0 = tic;
    roi = 256; mapT = 1.0; [~, iMap] = min(abs(tex/1e6 - mapT)); sigMap = [];
    for k = 1:n
        io = find(isOpen & tus == tex(k), 1); ic = find(~isOpen & tus == tex(k), 1);
        if isempty(io) || isempty(ic), continue; end
        [mO, vO] = cubeMedian(c.zipLin, names{io}, tmp, roi);
        [mC, vC] = cubeMedian(c.zipLin, names{ic}, tmp, roi);
        d = mO - mC;
        x = d(:); keep = true(size(x));
        for it = 1:3
            mu = mean(x(keep)); sd = std(x(keep)); keep = abs(x - mu) <= 5*sd;
        end
        S(k) = mean(x(keep)); Sstd(k) = std(x(keep)); Smed(k) = median(x);
        dark(k) = mean(mC(:)); satPix(k) = mean(mO(:) >= c.adcMax - 5);
        Vshot(k) = vO - vC;
        if k == iMap, sigMap = d; end
        fprintf('  %3d/%d  t = %.5f s  signal = %7.1f ADU  (%.0f s)\n', k, n, tex(k)/1e6, S(k), toc(t0));
    end
    t = tex / 1e6; ok = isfinite(S);
    sat = max(S(ok));
    % --- Fit 1: 5%-95% of saturation (Fig 4-19) ---
    r1 = ok & S >= 0.05*sat & S <= 0.95*sat;
    [p1, E1, res1] = linFit(t, S, r1);
    % --- Fit 2: < 2000 ADU (Fig 4-20) ---
    r2 = ok & S > 0 & S < 2000 & t > 0.001;
    [p2, E2, res2] = linFit(t, S, r2);
    % --- Fit 3: 2000-4000 ADU, no saturated points (Fig 4-21) ---
    r3 = ok & S >= 2000 & S < 0.98*sat;
    [p3, E3, res3] = linFit(t, S, r3);
    % --- Extra: photon transfer curve in central 512x512 ROI ---
    rp = ok & S > 50 & S < 0.85*sat;
    pp = polyfit(S(rp), Vshot(rp), 1); K_eADU = 1/pp(1);
    gainPTC = 1/(K_eADU*c.epsSi);               % ADU/eV equivalent
    fprintf('Saturation %.0f ADU; 5-95%% fit slope %.2f ADU/s, linearity error %.2f%%\n', sat, p1(1), E1);
    fprintf('<2000 ADU: slope %.2f ADU/s, E = %.2f%%;  2000-4000 ADU: slope %.2f ADU/s, E = %.2f%%\n', p2(1), E2, p3(1), E3);
    fprintf('Photon transfer: %.4f e-/ADU  (= %.4f ADU/eV)\n', K_eADU, gainPTC);
    OL = struct('t', t, 'S', S, 'Sstd', Sstd, 'Smed', Smed, 'dark', dark, 'satPix', satPix, 'Vshot', Vshot, ...
        'sat', sat, 'p1', p1, 'E1', E1, 'r1', r1, 'res1', res1, 'p2', p2, 'E2', E2, 'r2', r2, 'res2', res2, ...
        'p3', p3, 'E3', E3, 'r3', r3, 'res3', res3, 'K_eADU', K_eADU, 'gainPTC', gainPTC, 'ptc', pp, 'mapT', t(iMap));

    E = struct('C', c.E.C*gain, 'F', c.E.F*gain, 'MnKa', c.E.MnKa*gain, 'MnKb', c.E.MnKb*gain);
    % ---------------- Fig 4-19 style ----------------
    f = figure('Name', 'Optical linearity 5-95 pct', 'Color', 'w', 'Position', [40 40 1400 880]);
    tl = tiledlayout(f, 2, 1, 'TileSpacing', 'compact');
    nexttile;
    plot(t(ok), S(ok), 'o', 'Color', c.col.n1, 'MarkerSize', 5); hold on; grid on; box on;
    tt = linspace(0, max(t), 300); plot(tt, polyval(p1, tt), '--', 'Color', c.col.fit, 'LineWidth', 1.4);
    yline(sat, '-', 'Saturation', 'Color', c.col.orange, 'LineWidth', 1.4);
    t5 = (0.05*sat - p1(2))/p1(1); t95 = (0.95*sat - p1(2))/p1(1);
    xline([t5 t95], '--', {'5% Sat.', '95% Sat.'}, 'Color', c.col.spat, 'LineWidth', 1.3);
    xline(1.0, ':', 'Mode 1 | Mode 2', 'Color', [.4 .4 .4], 'LabelVerticalAlignment', 'bottom');
    ylabel('Mean Dark Subtracted Signal [ADU]'); xlabel('Exposure Time [Seconds]');
    title('IMX662: Linearity Measurements');
    legend({'Mean Signal', sprintf('Linear Fit Between 5%% and 95%% sat.: %.2f*x + %.2f', p1(1), p1(2))}, 'Location', 'southeast');
    nexttile;
    okr = ok & t > 0.001;   % the 16 us point has ~0 signal, so a percent residual is undefined
    plot(t(okr), 100*res1(okr)./polyval(p1, t(okr)), 'o', 'Color', c.col.n1, 'MarkerSize', 5); hold on; grid on; box on;
    xline([t5 t95], '--', {'5% Sat.', '95% Sat.'}, 'Color', c.col.spat); yline(0, 'k-');
    ylabel('Residuals [%]'); xlabel('Exposure Time [Seconds]');
    title(sprintf('IMX662: Linearity Residuals Between 5%% and 95%% Saturation   (linearity error, Eq. 3.4 = %.2f%%)', E1));
    legend({'Residuals'}, 'Location', 'southwest');
    imx662_caption(f, sprintf(['Optical linearity: the sensor was lit by a steady, uniform lamp (integrating sphere) and the exposure time was stepped from 2 s down to 0.02 s. ' ...
        'For each exposure the dark frame (shutter closed) is subtracted from the lit frame (shutter open). A perfectly linear sensor gives a straight line: double the exposure, double the signal. ' ...
        'Top: signal vs exposure with a straight-line fit using only points between 5%% and 95%% of saturation (green lines). Above ~%.0f ADU the pixels are full (saturated, orange line) and stop increasing. ' ...
        'Bottom: how far each point is from the line, in percent. The linearity error over this range is %.2f%%; the large values at the shortest exposures drive it (see next two plots).'], sat, E1), 0.17);

    % ---------------- Fig 4-20 / 4-21 style ----------------
    pieces = {r2, p2, E2, res2, 'below 2000 ADU', 'IMX662: Linearity Measurements Below 2000 ADU'; ...
              r3, p3, E3, res3, 'between 2000 and 4000 ADU', 'IMX662: Linearity Measurements Between 2000 and 4000 ADU'};
    for q = 1:2
        rr = pieces{q,1}; pq = pieces{q,2}; Eq = pieces{q,3}; rq = pieces{q,4};
        f = figure('Name', ['Optical linearity ' pieces{q,5}], 'Color', 'w', 'Position', [40 40 1400 880]);
        tl = tiledlayout(f, 2, 1, 'TileSpacing', 'compact');
        nexttile;
        plot(t(rr), S(rr), 'o', 'Color', c.col.n1, 'MarkerSize', 5); hold on; grid on; box on;
        tt = linspace(min(t(rr)), max(t(rr)), 200); plot(tt, polyval(pq, tt), '--', 'Color', c.col.fit, 'LineWidth', 1.4);
        if q == 1
            yline([E.C E.F], ':', {'C K\alpha (277 eV)', 'F K\alpha (677 eV)'}, 'Color', c.col.hot);
        else
            yline([E.MnKa E.MnKb], ':', {'Mn K\alpha (5900 eV)', 'Mn K\beta (6490 eV)'}, 'Color', c.col.hot);
        end
        ylabel('Mean Dark Subtracted Signal [ADU]'); xlabel('Exposure Time [Seconds]'); title(pieces{q,6});
        legend({'Mean Signal', sprintf('Linear Fit: %.2f*x + %.2f', pq(1), pq(2))}, 'Location', 'southeast');
        nexttile;
        plot(t(rr), 100*rq(rr)./polyval(pq, t(rr)), 'o', 'Color', c.col.n1, 'MarkerSize', 5); grid on; box on; yline(0, 'k-');
        ylabel('Residuals [%]'); xlabel('Exposure Time [Seconds]');
        title(sprintf('Residuals from the linear fit   (linearity error, Eq. 3.4 = %.2f%%)', Eq));
        imx662_caption(f, sprintf(['Same data, fitted separately for signals %s. Splitting the range lets each part be calibrated with its own straight line. ' ...
            'Dotted purple lines show where the lab''s X-ray lines fall in ADU (energy x gain), i.e. which linear regime matters for them. ' ...
            'Linearity error = (largest minus smallest residual) / (range of the fitted line) x 100 = %.2f%%: below ~0.5%% means the sensor is very linear in this range.'], pieces{q,5}, Eq), 0.17);
    end

    % ---------------- Extra: photon transfer curve ----------------
    f = figure('Name', 'Photon transfer curve', 'Color', 'w', 'Position', [40 40 1400 620]);
    tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact');
    nexttile; okp = ok & S > 0 & Vshot > 0;
    loglog(S(okp), Vshot(okp), 'o', 'Color', c.col.n1); hold on; grid on;
    ss = logspace(0, log10(sat), 100); loglog(ss, ss / K_eADU, '-', 'Color', c.col.fit, 'LineWidth', 1.4);
    xlabel('Mean Dark Subtracted Signal [ADU]'); ylabel('Shot-noise variance [ADU^2]');
    title(sprintf('Photon transfer curve: K = %.3f e^-/ADU', K_eADU));
    legend({'data (central 512\times512 pixels)', 'variance = signal / K (slope 1)'}, 'Location', 'northwest');
    nexttile; plot(S(ok), Vshot(ok), 'o', 'Color', c.col.n1); hold on; grid on;
    plot(ss, polyval(pp, ss), '-', 'Color', c.col.fit, 'LineWidth', 1.4);
    xlabel('Mean Dark Subtracted Signal [ADU]'); ylabel('Shot-noise variance [ADU^2]'); title('Same, linear axes (variance collapses at saturation)');
    imx662_caption(f, sprintf(['Extra check of the conversion gain using light instead of X-rays. Light arrives in random packets (photons), so the frame-to-frame scatter (variance) of a pixel grows in ' ...
        'proportion to its signal: variance = signal / K. The fitted K = %.3f electrons per ADU corresponds to %.4f ADU/eV, to compare with the Fe-55 X-ray gain. ' ...
        'The variance drops at the top because full pixels cannot fluctuate.'], K_eADU, gainPTC), 0.18);

    % ---------------- Extra: illumination map ----------------
    if ~isempty(sigMap)
        f = figure('Name', 'Optical illumination map', 'Color', 'w', 'Position', [40 40 1300 760]);
        imagesc(sigMap); axis image; colormap(f, imx662_viridis(256)); cb = colorbar; cb.Label.String = 'Dark Subtracted Signal [ADU]';
        v = sort(sigMap(:)); caxis([v(round(0.005*end)) v(round(0.995*end))]);
        xlabel('X Pixel'); ylabel('Y Pixel'); title(sprintf('IMX662 dark-subtracted signal at %.2f s exposure (median of 25 frames)', t(iMap)));
        set(gca, 'FontSize', c.fs);
        imx662_caption(f, ['Image of the lamp-lit sensor after dark subtraction. Shows how evenly the integrating sphere lit the sensor and any pixel-to-pixel response differences (PRNU). ' ...
            'A smooth gradient is illumination; speckle is pixel response non-uniformity.'], 0.13);
    end
end

function [p, E, res] = linFit(t, S, r)
    p = polyfit(t(r), S(r), 1); res = S - polyval(p, t);
    l = polyval(p, t(r)); E = (max(res(r)) - min(res(r))) / (max(l) - min(l)) * 100;
end

function [m, vTemp] = cubeMedian(zipPath, entry, tmp, roi)
    f = imx662_zip_extract(zipPath, entry, tmp); cl = onCleanup(@() delete(f));
    cube = imx662_read_cube(f);
    m = median(cube, 3);
    [ny, nx, ~] = size(cube);
    r = round(ny/2) + (-roi+1:roi); q = round(nx/2) + (-roi+1:roi);
    d = diff(cube(r, q, :), 1, 3); d = d - mean(mean(d, 1), 2);
    vTemp = mean(var(reshape(d, [], size(d,3)), 0, 1)) / 2;
end
