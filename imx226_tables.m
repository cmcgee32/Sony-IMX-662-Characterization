function imx226_tables(OUT)
%IMX226_TABLES  IMX226 X-ray linearity, resolution budget and summary tables (CSV).
%   Carissma (2026) next to Lupo (2025) for every quantity the thesis reports.
    c = imx662_cfg();
    load(fullfile(OUT, 'imx226_po210.mat'), 'Rpo'); load(fullfile(OUT, 'imx226_fe55.mat'), 'Rfe');
    S = load(fullfile(OUT, 'imx226_readnoise.mat')); RN = S.RN;
    XL = imx662_xray_linearity(Rfe, Rpo, RN);
    save(fullfile(OUT, 'imx226_xray_linearity.mat'), 'XL');
    ME = 'Carissma (2026)'; LU = 'Lupo (2025)';
    lines = {Rpo, 'C'; Rpo, 'F'; Rfe, 'MnKa'; Rfe, 'MnKb'};
    nm = {'C K-alpha', 'F K-alpha', 'Mn K-alpha', 'Mn K-beta'}; Eexp = [277 677 5900 6490];
    lu1 = {'47.2 +- 0.8', '51.1 +- 0.3', '142.9 +- 0.2', '154.5 +- 0.6'};
    lu4 = {'52.4 +- 0.7', '57.9 +- 0.3', '167.2 +- 0.5', '185.1 +- 1.2'};
    luC1 = [4349800 16752000 227150000 28291000]; luC4 = [7657200 32977000 805980000 112690000];
    luR1 = [0.1244 0.4738 56.1643 7.1084]; luR4 = [0.2156 0.9364 200.2699 27.5173];
    T3 = {}; T4 = {}; T5 = {}; T7 = {};
    for i = 1:4
        R = lines{i,1}; q = R.L(strcmp({R.L.short}, lines{i,2}));
        T3(end+1,:) = {nm{i}, Eexp(i), q.n1.centerEV, sprintf('%.1f +- %.1f', q.n1.fwhmEV, q.n1.fwhmEVerr), q.fano, lu1{i}}; %#ok<AGROW>
        T4(end+1,:) = {nm{i}, Eexp(i), q.n1234.centerEV, sprintf('%.1f +- %.1f', q.n1234.fwhmEV, q.n1234.fwhmEVerr), q.fano, lu4{i}}; %#ok<AGROW>
        T5(end+1,:) = {nm{i}, q.n1.area, luC1(i), q.n1234.area, luC4(i)}; %#ok<AGROW>
        T7(end+1,:) = {nm{i}, q.rate1, luR1(i), q.rate1234, luR4(i)}; %#ok<AGROW>
    end
    hdrE = {'Emission line', 'Expected energy [eV]', 'Measured center [eV]', ['FWHM, ' ME ' [eV]'], 'Fano limit FWHM [eV]', [LU ' FWHM [eV]']};
    writecell([hdrE; T3], fullfile(OUT, 'Table_IMX226_energy_resolution_single_pixel.csv'));
    writecell([hdrE; T4], fullfile(OUT, 'Table_IMX226_energy_resolution_multi_pixel.csv'));
    writecell([{'Emission line', ['Total counts single-pixel, ' ME], [LU ' single'], ['Total counts multi-pixel, ' ME], [LU ' multi']}; T5], ...
              fullfile(OUT, 'Table_IMX226_total_counts.csv'));
    writecell([{'Emission line', ['Count rate single-pixel, ' ME ' [Counts/Sec/mm^2]'], [LU ' single'], ['Count rate multi-pixel, ' ME ' [Counts/Sec/mm^2]'], [LU ' multi']}; T7], ...
              fullfile(OUT, 'Table_IMX226_count_rates.csv'));
    writecell({'Source', ['Hot pixels, ' ME], ['Percent of sensor, ' ME ' [%]'], 'Criterion', LU; ...
               'Po-210', Rpo.nHot, Rpo.hotPct, 'lab per-event hot-pixel flags', '66708 (0.56%)'; ...
               'Fe-55', Rfe.nHot, Rfe.hotPct, '>= 11 events in Fe-55 data cube 1 (Poisson tail)', '66312 (0.55%)'; ...
               'Fe-55 (Gaussian + 3 sigma of counts per pixel)', Rfe.nHotGauss, 100*Rfe.nHotGauss/12e6, sprintf('> %.1f events', Rfe.hotThr), 'threshold 84 counts'}, ...
              fullfile(OUT, 'Table_IMX226_hot_pixels.csv'));
    writecell({'Parameter', [ME ' [e- RMS]'], [LU ' [e- RMS]']; 'Average', RN.avg, 2.47; 'RMS', RN.rms, 3.23; 'Mode', RN.mode, 2.28; 'Median', RN.median, NaN}, ...
              fullfile(OUT, 'Table_IMX226_readout_noise.csv'));
    rows = {'Conversion gain from F K-alpha [ADU/eV]', Rpo.gain, '0.3654';
            'Conversion gain from Mn K-alpha [ADU/eV]', Rfe.gainMn, '-';
            'Conversion gain from X-ray linearity slope [ADU/eV]', XL.fit(1), '0.37216';
            'X-ray linearity error [%]', XL.linErr, '0.512'};
    if isfield(Rfe.GV, 'pct')
        rows(end+1:end+3,:) = {'Pixel-to-pixel gain variation, Mn K-alpha [%]', Rfe.GV.pct, 'not measured';
                               'Gain-variation broadening sigma at 5900 eV [eV]', Rfe.GV.sigmaEV, 'not measured';
                               'Pixels used for gain variation', Rfe.GV.nPix, '-'};
    end
    if exist(fullfile(OUT, 'imx226_optical_linearity.mat'), 'file')
        S = load(fullfile(OUT, 'imx226_optical_linearity.mat')); OL = S.OL;
        rows(end+1:end+5,:) = {'Optical linearity error, 5-95% saturation [%]', OL.E1, '-';
                               'Optical linearity error, Mode 1 (< 1000 ADU) [%]', OL.E2, '0.22';
                               'Optical linearity error, Mode 2 (1000-4000 ADU) [%]', OL.E3, '0.47';
                               'Optical saturation level [ADU]', OL.sat, '~4000';
                               'Conversion gain from optical photon transfer [ADU/eV]', OL.gainPTC, '-'};
    end
    if exist(fullfile(OUT, 'imx226_dark.mat'), 'file')
        S = load(fullfile(OUT, 'imx226_dark.mat')); DC = S.DC;
        rows(end+1:end+5,:) = {'Dark current, lab method, full detector [e-/pixel/s]', DC.fits(1,1), '0.43';
                               'Dark current, lab method, center region [e-/pixel/s]', DC.fits(2,1), '0.22';
                               'Dark current, lab method, edge rows [e-/pixel/s]', DC.fits(3,1), '0.75';
                               'Dark current, shot-noise method [e-/pixel/s]', DC.DCshot, '-';
                               'Hot dark-current pixels', DC.nHot, '-'};
    end
    writecell([{'Quantity', ME, [LU ' IMX226']}; rows], fullfile(OUT, 'Table_IMX226_summary_gain_linearity_dark.csv'));
    B = XL.budget; luU1 = containers.Map({'C K-alpha','F K-alpha','Mn K-alpha','Mn K-beta'}, {14.7, 10.8, 33.7, 38.9});
    luU4 = containers.Map({'C K-alpha','F K-alpha','Mn K-alpha','Mn K-beta'}, {17.6, 15.9, 50.0, 58.2});
    lu = nan(height(B), 1);
    for i = 1:height(B)
        k = strrep(strrep(B.Line{i}, 'C K-alpha', 'C K-alpha'), ' ', ' ');
        if isKey(luU1, k), if strcmp(B.Events{i}, 'n1'), lu(i) = luU1(k); else, lu(i) = luU4(k); end, end
    end
    B.Lupo2025_unaccounted_eV = lu;
    writetable(B, fullfile(OUT, 'Table_IMX226_energy_resolution_budget.csv'));
    fprintf('IMX226 tables written to %s\n', OUT);
end
