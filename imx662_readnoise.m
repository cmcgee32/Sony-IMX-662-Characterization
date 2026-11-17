function RN = imx662_readnoise(gain)
%IMX662_READNOISE  Readout noise per pixel (lab method, Lupo 2025 Sec. 3.2).
%   Every X-ray data chunk has a *_biasmap.fits with a per-pixel NOISE map
%   (std of 100 bias frames, ADU). The per-pixel median over all maps from the
%   X-ray runs gives one noise map, converted to electrons:
%       e- = ADU / gain[ADU/eV] / 3.65 eV
%   Reports average, RMS and mode of the distribution (Table 4-9 style).
    c = imx662_cfg();
    if nargin < 1 || isempty(gain), gain = c.gainRef; end
    tmp = fullfile(tempdir, 'imx662_rn'); if exist(tmp,'dir'), rmdir(tmp,'s'); end; mkdir(tmp);
    cln = onCleanup(@() rmdir(tmp, 's'));
    zips = {c.zipFe, c.zipPo}; maxPer = 46;
    list = {};
    for z = 1:2
        nm = imx662_zip_list(zips{z}, '_biasmap\.fits$'); nm = nm(~contains(nm, 'OlderFiles'));
        ii = round(linspace(1, numel(nm), min(maxPer, numel(nm))));
        for i = ii, list(end+1, :) = {zips{z}, nm{i}}; end %#ok<AGROW>
    end
    n = size(list, 1);
    Nst = zeros(c.ny, c.nx, n, 'uint16'); Bst = Nst; medN = zeros(n,1); medB = medN;
    fprintf('\n=== Readout noise: %d noise maps ===\n', n); t0 = tic;
    for k = 1:n
        f = imx662_zip_extract(list{k,1}, list{k,2}, tmp);
        B = fitsread(f, 'image', 1); N = fitsread(f, 'image', 2); delete(f);
        Nst(:,:,k) = uint16(round(N*100)); Bst(:,:,k) = uint16(round(B*10));
        medN(k) = median(N(:)); medB(k) = median(B(:));
        if mod(k, 10) == 1 || k == n, fprintf('  map %d/%d (%.0f s)\n', k, n, toc(t0)); end
    end
    noiseADU = double(median(Nst, 3)) / 100;  clear Nst
    biasADU  = double(median(Bst, 3)) / 10;   clear Bst
    noiseE = noiseADU / gain / c.epsSi;
    v = noiseE(:);
    be = 0:0.005:12; hc = histcounts(v, be); bc = be(1:end-1) + 0.0025;
    [~, im] = max(hc);
    RN.avg = mean(v); RN.rms = sqrt(mean(v.^2)); RN.mode = bc(im);
    RN.median = median(v); RN.noiseE = noiseE; RN.noiseADU = noiseADU; RN.biasADU = biasADU;
    RN.medN = medN; RN.medB = medB; RN.src = list(:,2); RN.gain = gain;
    fprintf('Readout noise [e- RMS]: average %.3f, RMS %.3f, mode %.3f, median %.3f\n', RN.avg, RN.rms, RN.mode, RN.median);

    % ---- Figure: noise histogram (lab-standard) ----
    f = figure('Name', 'Readout noise histogram', 'Color', 'w', 'Position', [40 40 1300 650]);
    stairs(bc, max(hc, 0.5), 'Color', c.col.n1, 'LineWidth', 0.7); set(gca, 'YScale', 'log'); grid on; box on; hold on;
    xline(RN.mode, '--', sprintf('mode %.2f e^-', RN.mode), 'Color', c.col.fit, 'LineWidth', 1.3);
    xline(RN.avg, ':', sprintf('average %.2f e^-', RN.avg), 'Color', c.col.reb, 'LineWidth', 1.3, 'LabelVerticalAlignment', 'middle');
    xline(RN.rms, '-.', sprintf('RMS %.2f e^-', RN.rms), 'Color', c.col.hot, 'LineWidth', 1.3, 'LabelVerticalAlignment', 'bottom');
    xlim([0 10]); ylim([0.5 max(hc)*3]);
    xlabel('Electron rms'); ylabel('Counts (number of pixels)'); title('IMX662 Readout Noise');
    legend({'IMX662 - Readout Noise'}, 'Location', 'northeast'); set(gca, 'FontSize', c.fs);
    imx662_caption(f, sprintf(['Readout noise = the random electronic "fuzz" added every time a pixel is read, in electrons (e-). Before each data chunk the camera took 100 dark frames and ' ...
        'measured each pixel''s noise; the per-pixel median over %d chunks is histogrammed here (log y). Lower is better: it sets the faintest X-ray that can be seen. ' ...
        'Mode (most common value) = %.2f e-, average = %.2f e-, RMS = %.2f e-. The tail to the right is a small number of noisier pixels. ' ...
        'The saw-tooth ripple is the same pattern reported in Lupo (2025); it is believed to come from the camera''s analog-to-digital converter.'], n, RN.mode, RN.avg, RN.rms), 0.2);

    % ---- Extra: noise map and bias map ----
    f = figure('Name', 'Readout noise and bias maps', 'Color', 'w', 'Position', [40 40 1500 560]);
    tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact');
    nexttile; imagesc(noiseE); axis image; colormap(gca, imx662_viridis(256)); cb = colorbar; cb.Label.String = 'Readout noise [e^- rms]';
    caxis([prctileV(v, 1) prctileV(v, 99)]); xlabel('X pixel'); ylabel('Y pixel'); title('Per-pixel readout noise map');
    nexttile; imagesc(biasADU); axis image; colormap(gca, imx662_viridis(256)); cb = colorbar; cb.Label.String = 'Bias level [ADU]';
    bv = biasADU(:); caxis([prctileV(bv, 1) prctileV(bv, 99)]); xlabel('X pixel'); ylabel('Y pixel'); title('Per-pixel bias (offset) map');
    imx662_caption(f, ['Left: readout noise of every pixel (electrons). Right: the bias, i.e. the signal each pixel reports with no light at all (about 100 ADU, an offset the camera ' ...
        'adds on purpose and that is subtracted from every frame). Stripes or blocks would point to readout-circuit structure; a uniform map means no such problems.'], 0.15);

    % ---- Extra: stability over the X-ray runs ----
    f = figure('Name', 'Bias and noise stability', 'Color', 'w', 'Position', [40 40 1200 520]);
    yyaxis left; plot(1:n, medB, 'o-', 'Color', c.col.reb); ylabel('Median bias [ADU]');
    yyaxis right; plot(1:n, medN / gain / c.epsSi, 's-', 'Color', c.col.fit); ylabel('Median readout noise [e^- rms]');
    grid on; xlabel('Data chunk (first 46: Fe-55 run, next 46: Po-210 run, ~5 min each)');
    xline(sum(contains(list(:,1), 'Fe55')) + 0.5, '--k', 'Fe-55 | Po-210');
    title('Bias and readout noise over the X-ray runs'); set(gca, 'FontSize', c.fs);
    imx662_caption(f, 'Median bias and median readout noise of the whole sensor for each data chunk. Flat lines mean the electronics were stable during the long X-ray runs.', 0.15);
end

function v = prctileV(x, p)
    x = sort(x(:)); v = x(max(1, min(numel(x), round(p/100*numel(x)))));
end
