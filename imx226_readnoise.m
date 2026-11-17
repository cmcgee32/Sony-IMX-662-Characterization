function RN = imx226_readnoise(folder, gain)
%IMX226_READNOISE  IMX226 readout noise from the lab bias maps (Lupo 2025, Sec. 3.2).
%   The Po-210 run has 353 bias maps (100 dark frames each, BIAS + NOISE
%   images). Nine maps spread evenly over runs 2-353 were streamed from Box;
%   for every pixel the median NOISE over the nine maps is used (lab method:
%   per-pixel median), converted ADU -> e- with the measured gain and 3.65 eV.
%   folder: med.u16 (per-pixel median noise x 200 [ADU]), bias4.f32 (bias of
%   the first map in 4x4 clusters [ADU]), meta.json (runs used, per-map
%   median bias and noise).
    c = imx662_cfg();
    if nargin < 2 || isempty(gain), gain = 0.3654; end
    M = jsondecode(fileread(fullfile(folder, 'meta.json')));
    fid = fopen(fullfile(folder, 'med.u16'), 'r', 'ieee-le'); q = fread(fid, Inf, '*uint16'); fclose(fid);
    fid = fopen(fullfile(folder, 'bias4.f32'), 'r', 'ieee-le'); b4 = fread(fid, Inf, '*single'); fclose(fid);
    noiseADU = reshape(double(q) / 200, 4000, 3000)';            % rows = Y
    biasADU = reshape(double(b4), 1000, 750)';
    noiseE = noiseADU / gain / c.epsSi; v = noiseE(:);
    be = 0:0.005:12; hc = histcounts(v, be); bc = be(1:end-1) + 0.0025; [~, im] = max(hc);
    RN = struct('avg', mean(v), 'rms', sqrt(mean(v.^2)), 'mode', bc(im), 'median', median(v), 'gain', gain, ...
                'runs', M.runs, 'perMap', M.perMap, 'bc', bc, 'hc', hc);
    fprintf('IMX226 readout noise [e- RMS]: average %.3f, RMS %.3f, mode %.3f, median %.3f\n', RN.avg, RN.rms, RN.mode, RN.median);
    lupo = [2.47 3.23 2.28];

    f = figure('Name', 'IMX226 readout noise histogram', 'Position', [40 40 1300 650]);
    stairs(bc, max(hc, 0.5), 'Color', c.col.fit, 'LineWidth', 0.7); set(gca, 'YScale', 'log'); hold on;
    xline(RN.mode, '--', sprintf('Mode %.2f e^-', RN.mode), 'Color', c.col.all, 'LineWidth', 1.3);
    xline(RN.avg, ':', sprintf('Average %.2f e^-', RN.avg), 'Color', c.col.reb, 'LineWidth', 1.3, 'LabelVerticalAlignment', 'middle');
    xline(RN.rms, '-.', sprintf('RMS %.2f e^-', RN.rms), 'Color', c.col.hot, 'LineWidth', 1.3, 'LabelVerticalAlignment', 'bottom');
    xlim([0 10]); ylim([0.5 max(hc)*3]);
    xlabel('Electron rms'); ylabel('Counts (number of pixels)'); title('IMX226 Readout Noise');
    legend({'IMX226 - Readout Noise'}, 'Location', 'northeast');
    imx662_style(f);

    % comparison with the IMX662 (Carissma 2026) if available
    f662 = fullfile(c.OUT, 'readnoise_results.mat');
    if exist(f662, 'file')
        S = load(f662); if isfield(S, 'RN'), R6 = S.RN; else, R6 = S; end
        v6 = R6.noiseE(:); h6 = histcounts(v6, be);
        f = figure('Name', 'IMX226 vs IMX662 readout noise', 'Position', [40 40 1300 650]);
        stairs(bc, max(hc/sum(hc), 1e-9), 'Color', c.col.fit, 'LineWidth', 0.9); hold on;
        stairs(bc, max(h6/sum(h6), 1e-9), 'Color', c.col.n1, 'LineWidth', 0.9);
        set(gca, 'YScale', 'log'); xlim([0 10]); ylim([1e-7 1e-1]);
        xlabel('Electron rms'); ylabel('Fraction of Pixels'); title('Readout Noise: IMX226 and IMX662');
        legend({sprintf('IMX226 (mode %.2f e^-)', RN.mode), sprintf('IMX662 (mode %.2f e^-)', R6.mode)}, 'Location', 'northeast');
        imx662_style(f);
    end

    f = figure('Name', 'IMX226 readout noise and bias maps', 'Position', [40 40 1500 560]);
    tiledlayout(1, 2, 'TileSpacing', 'compact');
    nexttile; nc = squeeze(mean(mean(reshape(noiseE, 4, 750, 4, 1000), 1), 3));
    imagesc([0 3999], [0 2999], nc); axis image; colormap(gca, imx662_viridis(256)); cb = colorbar; cb.Label.String = 'Readout noise [e^- rms]';
    caxis(prctile(nc(:), [1 99])); xlabel('X Pixel'); ylabel('Y Pixel'); title('IMX226 Readout Noise Map (4\times4-pixel clusters)');
    nexttile; imagesc([0 3999], [0 2999], biasADU); axis image; colormap(gca, imx662_viridis(256)); cb = colorbar; cb.Label.String = 'Bias level [ADU]';
    caxis(prctile(biasADU(:), [1 99])); xlabel('X Pixel'); ylabel('Y Pixel'); title('IMX226 Bias Map (4\times4-pixel clusters)');
    imx662_style(f);

    pm = M.perMap;
    f = figure('Name', 'IMX226 bias and noise stability', 'Position', [40 40 1200 520]);
    yyaxis left; plot(pm(:,1), pm(:,2), 'o-', 'Color', c.col.reb); ylabel('Median Bias [ADU]');
    yyaxis right; plot(pm(:,1), pm(:,3) / gain / c.epsSi, 's-', 'Color', c.col.fit); ylabel('Median Readout Noise [e^- rms]');
    xlabel('Po-210 Run Number (5 min each)'); title('IMX226 Bias and Readout Noise over the Po-210 Run');
    imx662_style(f);
    RN.lupo = lupo;
end
