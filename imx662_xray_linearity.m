function XL = imx662_xray_linearity(Fe, Po, RN)
%IMX662_XRAY_LINEARITY  X-ray linearity (Lupo 2025 Sec. 3.4.1) and energy-resolution budget (Sec. 3.6).
%   Peak centers [ADU] of C, F (Po-210) and Si K-alpha, Si escape, Mn K-alpha,
%   Mn K-beta (Fe-55), single-pixel fits, vs expected energy [eV]. Weighted
%   linear fit (weights 1/e^2, e = sigma/sqrt(area), Eq. 3.3). Residuals in ADU.
%   Resolution budget: R = 2.3548*sqrt(sigFano^2 + N*sigRead^2 + sigGain^2) (Eq. 3.5).
    c = imx662_cfg();
    sn = 'IMX662'; if isfield(Fe, 'sensor'), sn = Fe.sensor; end
    P = [];
    srcs = {Po, Fe};
    for s = 1:2
        L = srcs{s}.L;
        for j = 1:numel(L)
            q = L(j).n1;
            if ~isfinite(q.mu) || q.area < 50 || strcmp(L(j).short, 'O'), continue; end   % lab set: C, F, Si Ka, Si esc, Mn Ka, Mn Kb
            P(end+1).name = L(j).name; P(end).short = L(j).short; P(end).E = L(j).E; %#ok<AGROW>
            P(end).adu = q.mu; P(end).e = q.e33; P(end).area = q.area;
        end
    end
    E = [P.E]'; A = [P.adu]'; w = 1 ./ max([P.e]', 1e-6).^2;
    X = [E ones(size(E))]; W = diag(w); b = (X'*W*X) \ (X'*W*A);
    res = A - X*b; lfit = X*b;
    linErr = (max(res) - min(res)) / (max(lfit) - min(lfit)) * 100;
    fprintf('\nX-ray linearity: gain (slope) %.5f ADU/eV, offset %.2f ADU, linearity error %.3f%%\n', b(1), b(2), linErr);
    XL = struct('P', P, 'fit', b, 'res', res, 'linErr', linErr);

    f = figure('Name', 'X-ray linearity', 'Color', 'w', 'Position', [40 40 1300 860]);
    tl = tiledlayout(f, 2, 1, 'TileSpacing', 'compact');
    nexttile;
    plot(E, A, 'o', 'Color', c.col.n1, 'MarkerSize', 8, 'LineWidth', 1.5); hold on; grid on; box on;
    xx = linspace(0, 7000, 100); plot(xx, b(1)*xx + b(2), ':', 'Color', c.col.orange, 'LineWidth', 1.6);
    for i = 1:numel(P), text(E(i) + 80, A(i), P(i).name, 'FontSize', 9); end
    xlabel('Expected Energy Peak [eV]'); ylabel('Measured Energy Peak [ADU]'); title([sn ': X-Ray Linearity']);
    legend({'Energy Peak Positions', sprintf('Line of Best Fit: %.5f*x + %.2f', b(1), b(2))}, 'Location', 'northwest');
    nexttile;
    plot(E, res, 'o', 'Color', c.col.n1, 'MarkerSize', 8, 'LineWidth', 1.5); grid on; box on; yline(0, 'k-');
    xlabel('Expected Energy Peak [eV]'); ylabel('Residuals [ADU]');
    title(sprintf('%s: X-Ray Linearity Residuals (linearity error = %.3f%%)', sn, linErr)); legend({'Residuals'});
    imx662_caption(f, sprintf(['X-ray linearity: each emission line has a known energy (x) and the sensor measured its peak position in ADU (y, center of the single-pixel Gaussian fit). ' ...
        'If the sensor is linear the points lie on a straight line; the slope is the conversion gain (%.5f ADU/eV, vs %.5f from Mn K-alpha alone). ' ...
        'Bottom: distance of each point from the line in ADU. Linearity error = %.3f%%; values well below 1%% mean energies can be converted with one straight line from 277 eV to 6490 eV.'], ...
        b(1), Fe.gain, linErr), 0.18);

    % ---------------- Energy-resolution budget ----------------
    sigRead = RN.mode * c.epsSi;                     % eV, mode of readout noise
    rows = {};
    for s = 1:2
        L = srcs{s}.L;
        for j = 1:numel(L)
            if strcmp(L(j).short, 'O') || strcmp(L(j).short, 'SiEsc') || strcmp(L(j).short, 'SiKa'), continue; end
            sF = 2.3548*sqrt(c.fano*c.epsSi*L(j).E) / 2.3548;
            for cls = {'n1', 'n1234'}
                q = L(j).(cls{1});
                if strcmp(cls{1}, 'n1'), Np = 1; else, Np = max(L(j).meanNpix, 1); end
                Rcalc = 2.3548*sqrt(sF^2 + Np*sigRead^2);
                sG = 0; if strcmp(L(j).short, 'MnKa') && isfield(Fe.GV, 'sigmaEV'), sG = Fe.GV.sigmaEV; end
                RcalcG = 2.3548*sqrt(sF^2 + Np*sigRead^2 + sG^2);
                unacc = sqrt(max(q.fwhmEV^2 - Rcalc^2, 0)) / 2.3548;
                rows(end+1, :) = {strrep(L(j).name, '\alpha', '-alpha'), cls{1}, q.fwhmEV, q.fwhmEVerr, 2.3548*sF, Np, Rcalc, RcalcG, unacc}; %#ok<AGROW>
                rows{end,1} = strrep(rows{end,1}, '\beta', '-beta');
            end
        end
    end
    XL.budget = cell2table(rows, 'VariableNames', {'Line', 'Events', 'Measured_FWHM_eV', 'Measured_err_eV', 'Fano_FWHM_eV', ...
        'Mean_pixels', 'Calc_FWHM_Fano_Read_eV', 'Calc_FWHM_Fano_Read_Gain_eV', 'Unaccounted_sigma_eV'});
    XL.sigReadEV = sigRead;
    disp(XL.budget);

    f = figure('Name', 'Energy resolution budget', 'Color', 'w', 'Position', [40 40 1300 650]);
    m1 = strcmp(rows(:,2), 'n1');
    labs = rows(m1, 1); meas = cell2mat(rows(m1, 3)); fano = cell2mat(rows(m1, 5)); calc = cell2mat(rows(m1, 7)); calcG = cell2mat(rows(m1, 8));
    bh = bar(categorical(labs, labs), [fano calc calcG meas], 'grouped');
    bh(1).FaceColor = [0.75 0.75 0.75]; bh(2).FaceColor = c.col.n1; bh(3).FaceColor = c.col.hot; bh(4).FaceColor = c.col.orange;
    grid on; ylabel('Energy Resolution FWHM [eV]');
    legend({'Fano Limit', '+ Readout Noise', '+ Readout Noise + Gain Variation', 'Measured'}, 'Location', 'northwest');
    title([sn ' Single-Pixel Energy Resolution: Measured and Calculated']); set(gca, 'FontSize', c.fs);
    imx662_caption(f, sprintf(['Energy resolution (peak width, FWHM, smaller = better) for each line. Grey: the physical limit set by the statistics of electron creation in silicon (Fano limit). ' ...
        'Blue: adding the readout noise (mode = %.2f e- = %.1f eV). Purple: also adding pixel-to-pixel gain variation (only measurable for Mn K-alpha). ' ...
        'Orange: what was measured. The gap between orange and the calculated bars is noise not explained by these sources (listed in the table as "Unaccounted sigma").'], RN.mode, sigRead), 0.17);
end
