function E = imx290_ecc_analysis()
%IMX290_ECC_ANALYSIS  Per-pixel gain calibration of the Sony IMX290 (Fe-55).
%   Re-analysis of the IMX290 "ecc" data set on Box (psyduck/ecc, 2023-12-24
%   Fe-55 long run, 10 fps, B. Schneider). Each chunk_*.fits file holds the
%   pulse-height (ALLSPPH) spectra of 625 pixels, 4-ADU bins, 0-2200 ADU;
%   16 chunks = 10,000 pixels. Folders *_<N>evts hold the same pixels built
%   from only N events per pixel.
%
%   Method (Energy Calibration by Correlation, as in Schneider et al. 2024,
%   JATIS 10(3) 038001): for every pixel the gain g [ADU/eV] is the value
%   that maximises the correlation between the pixel's spectrum and the
%   reference Fe-55 spectrum (src/*gcorr_bis.txt) placed at PH = g*E.
%   Outputs: figures (PNG + FIG), CSV tables and imx290_results.mat in
%   Downloads/IMX290_results. Lab units: ADU, eV, ADU/eV.
    c = imx662_cfg();
    D = fullfile(c.D, 'Box_data', 'ecc');
    OUT = fullfile(c.D, 'IMX290_results'); if ~exist(OUT, 'dir'), mkdir(OUT); end
    ben.gain = 0.294; ben.rn = 1.9;                       % Schneider et al. 2024
    ben.fw1 = [136.1 123.6]; ben.fwAll = [144.5 134.6];   % Mn K-alpha FWHM before/after [eV]
    fanoMn = 2.3548*sqrt(c.fano*c.epsSi*c.E.MnKa);

    % ---------------- read all pixels ----------------
    [pc, S, XY, hd] = readSet(fullfile(D, 'data', '20231224_spec4calib_highCounts'));
    nPix = size(S, 2); nEv = sum(S, 1);
    syn = readmatrix(fullfile(D, 'src', '20231224_084701_IMX290_Fe55_longRun_10fps_886_grades_gcorr_bis.txt'), ...
                     'FileType', 'text', 'NumHeaderLines', 1);
    Es = syn(:,1)*1000; Is = syn(:,2);                    % eV, relative intensity
    fprintf('IMX290 ECC: %d pixels, median %d events/pixel, %s fps, EXPTIME %s\n', nPix, round(median(nEv)), hd.FPS, hd.EXPTIME);

    % ---------------- ECC gains (all events) ----------------
    grid = 0.2850:0.00005:0.3020;                         % ADU/eV
    [g, rho] = eccGain(pc, S, Es, Is, grid);
    gPk = peakGain(pc, S, c.E.MnKa);                      % independent check: Mn K-alpha peak per pixel
    G = median(g);
    fprintf('ECC gain: median %.5f ADU/eV, spread %.3f %% (rms)\n', G, 100*std(g)/mean(g));

    % ---------------- spectra before / after correction ----------------
    eb = 0:4:8000; ec = eb(1:end-1) + 2;
    % single gain: every pixel's 4-ADU bins map to the same energies (13.6 eV bins), so keep native binning
    erw = (pc(:) / G)'; raw = sum(S, 2)';
    cor = hEnergy(pc, S, g, eb);
    Fr = fitLine(erw, raw, c.E.MnKa); Fc = fitLine(ec, cor, c.E.MnKa);
    Fr2 = fitLine(erw, raw, c.E.MnKb); Fc2 = fitLine(ec, cor, c.E.MnKb);
    fprintf('Mn Ka FWHM: %.1f -> %.1f eV;  Mn Kb: %.1f -> %.1f eV\n', Fr.fw, Fc.fw, Fr2.fw, Fc2.fw);

    % ---------------- calibration accuracy vs events per pixel ----------------
    Nlist = [10 50 100 200 500 600 800 1000 1200 1400];
    acc = nan(numel(Nlist), 5);   % N, rms dev [%], rms dev [eV @ 5.9 keV], FWHM corr [eV], FWHM err
    for k = 1:numel(Nlist)
        [~, Sk] = readSet(fullfile(D, 'data', sprintf('20231224_spec4calib_highCounts_%devts', Nlist(k))));
        gk = eccGain(pc, Sk, Es, Is, grid);
        d = (gk - g) ./ g;
        Fk = fitLine(ec, hEnergy(pc, S, gk, eb), c.E.MnKa);
        acc(k, :) = [Nlist(k), 100*rmsRobust(d), rmsRobust(d)*c.E.MnKa, Fk.fw, Fk.fwe];
        fprintf('  N = %4d: gain rms error %.3f %%, corrected FWHM %.1f eV\n', Nlist(k), acc(k,2), acc(k,4));
    end

    % ================= FIGURES (imx290_plot.m) =================
    E0 = struct('g', g, 'gPk', gPk, 'rho', rho, 'XY', XY, 'G', G, 'Fr', Fr, 'Fc', Fc, 'Fr2', Fr2, 'Fc2', Fc2, ...
                'acc', acc, 'ec', ec, 'erw', erw, 'raw', raw, 'cor', cor, 'ben', ben, 'nAll', median(nEv));
    pg = imx290_plot(E0);
    imx662_save_figs('IMX290', OUT);

    % ================= TABLES =================
    T1 = {'Conversion gain, median ECC [ADU/eV]', G, ben.gain;
          'Pixel-to-pixel gain variation (Gaussian sigma) [%]', 100*pg(3)/pg(2), NaN;
          'Gain variation, rms of all pixels [%]', 100*std(g)/mean(g), NaN;
          'Mn K-alpha FWHM, single gain [eV]', Fr.fw, ben.fw1(1);
          'Mn K-alpha FWHM, per-pixel gain [eV]', Fc.fw, ben.fw1(2);
          'Mn K-beta FWHM, single gain [eV]', Fr2.fw, NaN;
          'Mn K-beta FWHM, per-pixel gain [eV]', Fc2.fw, NaN;
          'Improvement in Mn K-alpha FWHM [eV]', Fr.fw - Fc.fw, ben.fw1(1) - ben.fw1(2);
          'Fano limit at Mn K-alpha [eV]', fanoMn, fanoMn;
          'Pixels calibrated', nPix, NaN;
          'Median events per pixel', median(nEv), NaN};
    writecell([{'Quantity', 'Carissma (2026)', 'Schneider et al. (2024)'}; T1], fullfile(OUT, 'Table_IMX290_ECC_summary.csv'));
    writecell([{'Events per pixel', 'RMS gain error [%]', 'RMS gain error at 5.9 keV [eV]', 'Corrected Mn K-alpha FWHM [eV]', 'FWHM error [eV]'}; num2cell(acc)], ...
              fullfile(OUT, 'Table_IMX290_accuracy_vs_events.csv'));
    E = struct('g', g, 'gPk', gPk, 'rho', rho, 'XY', XY, 'G', G, 'pg', pg, 'Fr', Fr, 'Fc', Fc, 'Fr2', Fr2, 'Fc2', Fc2, ...
               'acc', acc, 'ec', ec, 'erw', erw, 'raw', raw, 'cor', cor, 'ben', ben, 'nAll', median(nEv));
    save(fullfile(OUT, 'imx290_results.mat'), '-struct', 'E');
end

% ---------------------------------------------------------------- helpers
function [pc, S, XY, hd] = readSet(folder)
    fs = dir(fullfile(folder, 'chunk_*.fits'));
    k = cellfun(@(s) sscanf(s, 'chunk_%d'), {fs.name}); [~, o] = sort(k); fs = fs(o);
    S = []; XY = [];
    for i = 1:numel(fs)
        fn = fullfile(folder, fs(i).name);
        info = fitsinfo(fn); kw = info.BinaryTable(1).Keywords;
        T = fitsread(fn, 'binarytable');
        pc = double(T{2});
        S = [S, double(cell2mat(T(4:end)))]; %#ok<AGROW>
        nch = numel(T) - 3; xy = zeros(nch, 2);
        for j = 1:nch
            v = kw{strcmp(kw(:,1), sprintf('IDCH%03d', j-1)), 2};
            xy(j, :) = sscanf(strrep(v, ',', ' '), '%d %d')';
        end
        XY = [XY; xy]; %#ok<AGROW>
    end
    hd.FPS = num2str(kw{strcmp(kw(:,1), 'FPS'), 2}); hd.EXPTIME = num2str(kw{strcmp(kw(:,1), 'EXPTIME'), 2});
end

function [g, rhoMax] = eccGain(pc, S, Es, Is, grid)
    sel = pc > 1000 & pc < 2100;
    M = zeros(numel(grid), nnz(sel));
    for k = 1:numel(grid), M(k, :) = interp1(Es, Is, pc(sel)/grid(k), 'linear', 0); end
    M = M - mean(M, 2); M = M ./ vecnorm(M, 2, 2);
    O = S(sel, :); O = O - mean(O, 1); O = O ./ max(vecnorm(O, 2, 1), eps);
    C = M * O; [~, k] = max(C, [], 1); k = min(max(k, 2), numel(grid)-1);
    idx = sub2ind(size(C), k, 1:size(C,2));
    y0 = C(idx-1); y1 = C(idx); y2 = C(idx+1);
    d = 0.5*(y0 - y2) ./ (y0 - 2*y1 + y2); d(~isfinite(d)) = 0;
    g = grid(k) + d*(grid(2) - grid(1)); rhoMax = y1;
end

function g = peakGain(pc, S, E0)
    g = zeros(1, size(S, 2));
    for j = 1:size(S, 2)
        s = S(:, j); m = pc > 1550 & pc < 1950; [~, i] = max(s .* m); mu = pc(i);
        for it = 1:4
            w = abs(pc - mu) < 45; mu = sum(pc(w).*s(w)) / max(sum(s(w)), 1);
        end
        g(j) = mu / E0;
    end
end

function h = hEnergy(pc, S, g, eb)
    E = pc(:) ./ g(:)';                       % bins x pixels, eV
    W = S(:);  [~, ~, bi] = histcounts(E(:), eb); ok = bi > 0;
    h = accumarray(bi(ok), W(ok), [numel(eb)-1 1])';
end

function F = fitLine(ec, y, E0)
    m = ec > E0 - 250 & ec < E0 + 250; [~, i] = max(y .* m); mu = ec(i);
    m = ec > mu - 150 & ec < mu + 150;
    [p, pe] = imx662_gfit(ec(m), y(m), [max(y(m)) mu 55 min(y(m))], true);
    F.p = p; F.mu = p(2); F.fw = 2.3548*p(3); F.fwe = 2.3548*pe(3);
end

function r = rmsRobust(d)
    r = 1.4826 * median(abs(d - median(d)));
end
