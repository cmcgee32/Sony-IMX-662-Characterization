function BL = imx226_beamline(folder, gain)
%IMX226_BEAMLINE  IMX226 background data from the MIT polarimetry beamline (NE83, Nov 2024).
%   Box: psyduck/polarimetryBeamline/20241113, eight 300-frame cubes
%   "20241113_113533_IMX226_Bk_5fps_<10..17>_datacube.fits.gz" (5 fps, 0.2 s
%   frames, 23 C, analog gain 3). "Bk" = background: the detector on the
%   beamline with no X-ray line visible above threshold. The first 60 frames of
%   each 7.2 GB cube were streamed from Box and reduced in the browser:
%   bias = mean of frames 1-10; frames 11-60 searched for events with the
%   lab thresholds (event 40 ADU on a local-maximum pixel, split 9 ADU, 3x3).
%   Hot (warm) pixels: pixels that trigger in >= 2 frames of a cube.
%   Results: background event rate [Counts/Sec/cm^2], background spectrum,
%   spatial distribution, event sizes, bias stability between cubes, and a
%   comparison with the lab Po-210 bias / noise of the same sensor (2025).
    c = imx662_cfg();
    if nargin < 2 || isempty(gain), gain = 0.3654; end
    M = jsondecode(fileread(fullfile(folder, 'meta.json')));
    areaCM2 = 7.4*5.55/100; tf = 0.2;
    K = numel(M.cubes); all = []; rate = nan(K,1); nHotC = rate; biasMed = rate; nfr = rate; k0 = [];
    for i = 1:K
        q = M.cubes(i);
        fid = fopen(fullfile(folder, sprintf('cube_%s_events.i32', q.k)), 'r', 'ieee-le'); e = fread(fid, Inf, '*int32'); fclose(fid);
        e = double(reshape(e, 6, [])');
        pix = e(:,3)*4000 + e(:,2); [u, ~, j] = unique(pix); nper = accumarray(j, 1);
        hot = ismember(pix, u(nper >= 2)); nHotC(i) = nnz(nper >= 2);   % warm pixel: triggers in >= 2 of the 50 frames
        e = e(~hot, :); nf = q.nFrames - 10; nfr(i) = nf;   % frames searched for events (bias frames excluded)
        rate(i) = size(e, 1) / (nf*tf) / areaCM2;
        all = [all; e, repmat(str2double(q.k), size(e,1), 1)]; %#ok<AGROW>
        fid = fopen(fullfile(folder, sprintf('cube_%s_bias4.f32', q.k)), 'r', 'ieee-le'); b = fread(fid, Inf, '*single'); fclose(fid);
        biasMed(i) = median(double(b)); k0(i) = str2double(q.k); %#ok<AGROW>
        if i == 1, bias1 = reshape(double(b), 1000, 750)'; end
    end
    E = all(:,4) / gain;   % eV
    hiE = E > 300;   % above the noise-trigger peak: X-ray-like events
    rateHi = nan(K,1); for i = 1:K, rateHi(i) = nnz(hiE & all(:,7) == k0(i)) / (nfr(i)*tf) / areaCM2; end
    BL = struct('rate', rate, 'nHot', nHotC, 'biasMed', biasMed, 'cubes', k0, 'nEvents', size(all,1), ...
                'meanRate', size(all,1) / (sum(nfr)*tf) / areaCM2, 'gain', gain, 'tTot', sum(nfr)*tf, ...
                'rateHi', rateHi, 'meanRateHi', nnz(hiE) / (sum(nfr)*tf) / areaCM2, 'nHi', nnz(hiE));
    fprintf('IMX226 beamline background: %d cubes, %.0f s, %d events after hot-pixel removal -> %.3f Counts/Sec/cm^2\n', ...
            K, BL.tTot, BL.nEvents, BL.meanRate);
    fprintf('  events above 300 eV: %d -> %.3f Counts/Sec/cm^2\n', BL.nHi, BL.meanRateHi);

    f = figure('Name', 'IMX226 beamline background spectrum', 'Position', [40 40 1300 650]);
    eb = 0:10:3000; h1 = histcounts(E(all(:,5) == 1), eb); h4 = histcounts(E(all(:,5) <= 4), eb); h9 = histcounts(E, eb);
    ec = eb(1:end-1) + 5;
    stairs(ec, max(h9, 0.5), 'Color', c.col.all, 'LineWidth', 1.4); hold on;
    stairs(ec, max(h4, 0.5), 'Color', c.col.n1234, 'LineWidth', 1.2); stairs(ec, max(h1, 0.5), 'Color', c.col.n1, 'LineWidth', 1.2);
    set(gca, 'YScale', 'log'); xline(40/gain, '--', 'Event threshold (40 ADU)', 'Color', c.col.fit);
    xline([c.E.C c.E.O c.E.F], ':', {'C K\alpha', 'O K\alpha', 'F K\alpha'}, 'Color', [.5 .5 .5]);
    xlabel('Energy [eV]'); ylabel('Counts'); legend({'all', 'n1234', 'n1'}, 'Location', 'northeast');
    title(sprintf('IMX226 Polarimetry Beamline Background Spectrum (%d cubes, %.0f s)', K, BL.tTot));
    imx662_style(f);

    f = figure('Name', 'IMX226 beamline background map', 'Position', [40 40 1300 760]);
    cl = 100; H = histcounts2(all(:,3), all(:,2), 0:cl:3000, 0:cl:4000) / (BL.tTot * (cl*1.85e-4)^2);
    imagesc([0 3999], [0 2999], H); axis image; colormap(f, imx662_viridis(256)); cb = colorbar; cb.Label.String = 'Count Rate [Counts/Sec/cm^2]';
    xlabel('X Pixel'); ylabel('Y Pixel'); title('IMX226 Beamline Background Event Map (100\times100-pixel clusters)');
    imx662_style(f);

    f = figure('Name', 'IMX226 beamline background rate and bias', 'Position', [40 40 1300 520]);
    tiledlayout(1, 2, 'TileSpacing', 'compact');
    nexttile; plot(k0, rate, 'o-', 'Color', c.col.n1, 'MarkerFaceColor', c.col.n1); hold on;
    yline(BL.meanRate, '--', sprintf('Mean %.3f', BL.meanRate), 'Color', c.col.fit);
    plot(k0, rateHi, 's-', 'Color', c.col.n1234, 'MarkerFaceColor', c.col.n1234); set(gca, 'YScale', 'log');
    legend({'All events (> 40 ADU)', sprintf('Mean, all: %.1f', BL.meanRate), sprintf('Events > 300 eV (mean %.2f)', BL.meanRateHi)}, 'Location', 'east');
    xlabel('Beamline Cube Number'); ylabel('Background Rate [Counts/Sec/cm^2]'); title('Event Rate per Cube');
    nexttile; plot(k0, biasMed, 's-', 'Color', c.col.reb, 'MarkerFaceColor', c.col.reb); hold on;
    xlabel('Beamline Cube Number'); ylabel('Median Residual Offset [ADU]'); title('Residual Offset per Cube (camera bias already subtracted)');
    imx662_style(f);

    f = figure('Name', 'IMX226 beamline event size', 'Position', [40 40 900 600]);
    hn = histcounts(all(:,5), 0.5:1:9.5); bar(1:9, 100*hn/sum(hn), 'FaceColor', c.col.n1);
    xlabel('Number of Pixels in the Event'); ylabel('Percent of Events [%]'); title('IMX226 Beamline Background: Event Size Distribution');
    imx662_style(f);

    f = figure('Name', 'IMX226 beamline bias map', 'Position', [40 40 1300 760]);
    imagesc([0 3999], [0 2999], bias1); axis image; colormap(f, imx662_viridis(256)); cb = colorbar; cb.Label.String = 'Residual Offset [ADU]';
    caxis(prctile(bias1(:), [1 99])); xlabel('X Pixel'); ylabel('Y Pixel'); title('IMX226 Beamline Residual Offset Map, Cube 10 (4\times4-pixel clusters)');
    imx662_style(f);

    writecell([{'Cube', 'Frames used', 'Warm pixels (>= 2 triggers)', 'Event rate, all [Counts/Sec/cm^2]', 'Event rate above 300 eV [Counts/Sec/cm^2]', 'Median residual offset [ADU]'}; ...
               num2cell([k0(:) nfr nHotC rate rateHi biasMed])], fullfile(c.D, 'IMX226_results', 'Table_IMX226_beamline_background.csv'));
    BL.E = E; BL.events = all;
end
