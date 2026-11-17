function OL = imx226_optical_linearity(folder, gain)
%IMX226_OPTICAL_LINEARITY  IMX226 optical linearity (lab method, Lupo 2025 Sec. 3.4.2 / 4.4.2).
%   101 exposure pairs (shutter open / closed, 3.0 s down to 0.03 s) from
%   IMX226_Linearity_20250709.zip were streamed from Box (46 GB). For each
%   25-frame cube the first 8 frames were read; for each frame the 5-sigma
%   clipped mean over all pixels (3 iterations) was taken, and the median
%   over frames gives the cube level. Signal = open - closed (dark subtracted).
%   As in Lupo (2025), the IMX226 changes operating mode above 0.6 s, so the
%   fits are done separately for Mode 1 (t <= 0.6 s, below ~1000 ADU) and
%   Mode 2 (t > 0.6 s, 1000-4000 ADU, saturated points excluded), plus the
%   5-95 % of saturation fit. Linearity error: Eq. 3.4.
%   Photon transfer from the per-pixel temporal variance in a central
%   400 x 400-pixel region.
    c = imx662_cfg();
    if nargin < 2 || isempty(gain), gain = 0.3654; end
    M = jsondecode(fileread(fullfile(folder, 'meta.json')));
    C = M.cubes; tus = [C.tus]'; isOpen = logical([C.open]');
    lev = arrayfun(@(q) median(q.frames(:,1)), C); rv = [C.roiVar]'; rm = [C.roiMean]';
    tex = unique(tus); n = numel(tex); S = nan(n,1); dark = S; V = S; Sm = S;
    for k = 1:n
        io = find(isOpen & tus == tex(k), 1); ic = find(~isOpen & tus == tex(k), 1);
        if isempty(io) || isempty(ic), continue; end
        S(k) = lev(io) - lev(ic); dark(k) = lev(ic); V(k) = rv(io) - rv(ic); Sm(k) = rm(io) - rm(ic);
    end
    t = tex / 1e6; ok = isfinite(S);
    sat = max(S(ok));
    r1 = ok & S >= 0.05*sat & S <= 0.95*sat;
    r2 = ok & t < 0.6 - 1e-9 & S > 0;                  % Mode 1 (below the mode switch at 0.6 s, < ~1000 ADU)
    r3 = ok & t >= 0.6 - 1e-9 & S < 0.98*sat;          % Mode 2 (from 0.6 s, 1000-4000 ADU, not saturated)
    [p1, E1, res1] = linFit(t, S, r1); [p2, E2, res2] = linFit(t, S, r2); [p3, E3, res3] = linFit(t, S, r3);
    rp = ok & Sm > 50 & Sm < 0.8*sat & V > 0;
    pp = polyfit(Sm(rp), V(rp), 1); K = 1/pp(1); gPTC = 1/(K*c.epsSi);
    fprintf('IMX226 optical linearity: sat %.0f ADU; 5-95%%: E = %.2f%%; Mode 1 (<=0.6 s): E = %.2f%%; Mode 2: E = %.2f%%; PTC %.3f e-/ADU (%.4f ADU/eV)\n', ...
            sat, E1, E2, E3, K, gPTC);
    OL = struct('t', t, 'S', S, 'dark', dark, 'V', V, 'Sm', Sm, 'sat', sat, 'p1', p1, 'E1', E1, 'p2', p2, 'E2', E2, ...
                'p3', p3, 'E3', E3, 'K_eADU', K, 'gainPTC', gPTC, 'lupo', [0.22 0.47]);

    Eadu = struct('C', c.E.C*gain, 'F', c.E.F*gain);
    % ---- overall (Fig 4-22 style) ----
    f = figure('Name', 'IMX226 optical linearity', 'Position', [40 40 1400 880]);
    tiledlayout(2, 1, 'TileSpacing', 'compact');
    nexttile; plot(t(ok), S(ok), 'o', 'Color', c.col.n1, 'MarkerSize', 5); hold on;
    tt = linspace(0, max(t), 300); plot(tt, polyval(p1, tt), '--', 'Color', c.col.fit, 'LineWidth', 1.4);
    yline(sat, '-', 'Saturation', 'Color', c.col.orange, 'LineWidth', 1.4);
    xline(0.6, ':', 'Mode 1 | Mode 2', 'Color', c.col.fit, 'LineWidth', 1.4, 'LabelVerticalAlignment', 'bottom');
    ylabel('Mean Dark Subtracted Signal [ADU]'); xlabel('Exposure Time [Seconds]'); title('IMX226: Linearity Measurements');
    legend({'Mean Signal', sprintf('Linear Fit Between 5%% and 95%% sat.: %.2f*x + %.2f', p1(1), p1(2))}, 'Location', 'southeast');
    nexttile; okr = r1;   % residuals inside the 5-95 % fit range
    plot(t(okr), 100*res1(okr)./polyval(p1, t(okr)), 'o', 'Color', c.col.n1, 'MarkerSize', 5); hold on; yline(0, 'k-');
    xline(0.6, ':', 'Color', c.col.fit, 'LineWidth', 1.2);
    ylabel('Residuals [%]'); xlabel('Exposure Time [Seconds]'); title(sprintf('IMX226: Linearity Residuals (linearity error = %.2f%%)', E1));
    imx662_style(f);

    pieces = {r2, p2, E2, res2, 'IMX226: Linearity Measurements Below 1000 ADU (Mode 1)'; ...
              r3, p3, E3, res3, 'IMX226: Linearity Measurements Above 1000 ADU (Mode 2)'};
    for q = 1:2
        rr = pieces{q,1}; pq = pieces{q,2};
        f = figure('Name', sprintf('IMX226 optical linearity mode %d', q), 'Position', [40 40 1400 880]);
        tiledlayout(2, 1, 'TileSpacing', 'compact');
        nexttile; plot(t(rr), S(rr), 'o', 'Color', c.col.n1, 'MarkerSize', 5); hold on;
        tt = linspace(min(t(rr)), max(t(rr)), 200); plot(tt, polyval(pq, tt), '--', 'Color', c.col.fit, 'LineWidth', 1.4);
        if q == 1, yline([Eadu.C Eadu.F], ':', {'C K\alpha (277 eV)', 'F K\alpha (677 eV)'}, 'Color', c.col.hot); end
        ylabel('Mean Dark Subtracted Signal [ADU]'); xlabel('Exposure Time [Seconds]'); title(pieces{q,5});
        legend({'Mean Signal', sprintf('Linear Fit: %.2f*x + %.2f', pq(1), pq(2))}, 'Location', 'southeast');
        nexttile; plot(t(rr), 100*pieces{q,4}(rr)./polyval(pq, t(rr)), 'o', 'Color', c.col.n1, 'MarkerSize', 5); hold on; yline(0, 'k-');
        ylabel('Residuals [%]'); xlabel('Exposure Time [Seconds]');
        title(sprintf('Residuals from the Linear Fit (linearity error = %.2f%%)', pieces{q,3}));
        imx662_style(f);
    end

    f = figure('Name', 'IMX226 photon transfer curve', 'Position', [40 40 1400 620]);
    tiledlayout(1, 2, 'TileSpacing', 'compact');
    nexttile; okp = ok & Sm > 0 & V > 0; loglog(Sm(okp), V(okp), 'o', 'Color', c.col.n1); hold on;
    ss = logspace(0, log10(sat), 100); loglog(ss, ss/K, '-', 'Color', c.col.fit, 'LineWidth', 1.4);
    xlabel('Mean Dark Subtracted Signal [ADU]'); ylabel('Shot-Noise Variance [ADU^2]');
    title(sprintf('IMX226 Photon Transfer Curve: K = %.3f e^-/ADU', K));
    legend({'Central 400\times400 pixels', 'Variance = Signal / K'}, 'Location', 'northwest');
    nexttile; plot(Sm(ok), V(ok), 'o', 'Color', c.col.n1); hold on; plot(ss, polyval(pp, ss), '-', 'Color', c.col.fit, 'LineWidth', 1.4);
    xlabel('Mean Dark Subtracted Signal [ADU]'); ylabel('Shot-Noise Variance [ADU^2]'); title('Photon Transfer Curve (Linear Axes)');
    imx662_style(f);

    if isfield(M, 'coarse') && ~isempty(M.coarse)
        f = figure('Name', 'IMX226 optical illumination map', 'Position', [40 40 1300 760]);
        cm = reshape(M.coarse, 80, 60)';
        imagesc([0 3999], [0 2999], cm); axis image; colormap(f, imx662_viridis(256)); cb = colorbar; cb.Label.String = 'Signal [ADU]';
        xlabel('X Pixel'); ylabel('Y Pixel'); title(sprintf('IMX226 Illuminated Frame at %.2f s (50\\times50-pixel bins)', M.coarseT/1e6));
        imx662_style(f);
    end
end

function [p, E, res] = linFit(t, S, r)
    p = polyfit(t(r), S(r), 1); res = S - polyval(p, t);
    l = polyval(p, t(r)); E = (max(res(r)) - min(res(r))) / (max(l) - min(l)) * 100;
end
