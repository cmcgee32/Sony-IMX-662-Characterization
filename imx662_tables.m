function imx662_tables(Fe, Po, RN, XL, OL, DC)
%IMX662_TABLES  All result tables (equivalents of Lupo 2025 Tables 4-1 ... 4-15) as CSV + figures.
%   Each table has a "Lupo (2025)" column with the published IMX662 value for comparison.
    c = imx662_cfg(); out = c.OUT;
    nm = @(s) strrep(strrep(s, '\alpha', '-alpha'), '\beta', '-beta');
    getL = @(S, short) S.L(strcmp({S.L.short}, short));
    lines = {Po, 'C'; Po, 'F'; Fe, 'MnKa'; Fe, 'MnKb'};
    ref.res1 = [47.3 58.7 147.9 162.0]; ref.res1e = [0.4 0.3 1.5 1.8];
    ref.res4 = [49.4 61.1 151.2 165.5]; ref.res4e = [0.4 0.3 1.5 1.8]; ref.res4e = [0.6 0.3 1.5 1.8];
    ref.cnt1 = [539860 2492000 439180000 62186000]; ref.cnt4 = [562950 2607900 472530000 67076000];
    ref.rate1 = [0.0216 0.0995 164.8459 23.4055]; ref.rate4 = [0.0228 0.1047 177.7275 25.2484];
    ref.calc = [39.6 50.0 121.1 126.6]; ref.un1 = [11.0 13.1 36.1 42.9]; ref.un4 = [12.5 14.9 38.5 45.3];

    % ---------- Tables 4-1 / 4-2: hot pixels ----------
    D = {'Po-210', Po.nHot, Po.hotPct, sprintf('>= %g events per pixel', Po.hotThr), '44634 (2.15%)'; ...
         'Fe-55',  Fe.nHot, Fe.hotPct, sprintf('> %.1f events (center + 3 sigma)', Fe.hotThr), '15259 (0.74%)'};
    writecell([{'Source', 'Hot pixels', 'Percent of sensor [%]', 'Threshold', 'Lupo (2025)'}; D], fullfile(out, 'Table_01_hot_pixels.csv'));
    imx662_table_fig('Table hot pixels', 'Tables 4-1/4-2 equivalent: Hot pixels removed from the X-ray data', ...
        {'Radioactive source', 'Number of hot pixels', 'Percent of sensor [%]', 'Threshold used', 'Lupo (2025) IMX662'}, D, ...
        ['Hot pixels are pixels that report far more events than a real X-ray rate could explain (defects or noisy pixels); they are excluded from every spectrum, light curve and map. ' ...
         'Fe-55 uses the Gaussian center + 3 sigma rule on the events-per-pixel histogram, Po-210 uses >= 10 events (lab method). Percent = hot pixels / 2,073,600 pixels.']);

    % ---------- Tables 4-3 / 4-4: energy resolution ----------
    D1 = {}; D4 = {};
    for i = 1:4
        L = getL(lines{i,1}, lines{i,2});
        D1(i,:) = {nm(L.name), L.E, L.n1.centerEV, sprintf('%.1f +- %.1f', L.n1.fwhmEV, L.n1.fwhmEVerr), L.fano, sprintf('%.1f +- %.1f', ref.res1(i), ref.res1e(i))}; %#ok<AGROW>
        D4(i,:) = {nm(L.name), L.E, L.n1234.centerEV, sprintf('%.1f +- %.1f', L.n1234.fwhmEV, L.n1234.fwhmEVerr), L.fano, sprintf('%.1f +- %.1f', ref.res4(i), ref.res4e(i))}; %#ok<AGROW>
    end
    H = {'Emission line', 'Expected energy [eV]', 'Measured center [eV]', 'Energy resolution FWHM [eV]', 'Fano limit FWHM [eV]', 'Lupo (2025) FWHM [eV]'};
    writecell([H; D1], fullfile(out, 'Table_03_energy_resolution_single_pixel.csv'));
    writecell([H; D4], fullfile(out, 'Table_04_energy_resolution_multi_pixel.csv'));
    capR = ['Energy resolution = full width at half maximum (FWHM) of the Gaussian fit to each emission line, in eV: how well the sensor can tell two nearby X-ray energies apart (smaller is better). ' ...
            'The Fano limit is the best physically possible value for silicon (Eq. 3.1, F = 0.115, 3.65 eV per electron). Errors are 1-sigma fit uncertainties.'];
    imx662_table_fig('Table energy resolution single', 'Table 4-3 equivalent: Energy resolution for single-pixel events (n1)', H, D1, capR);
    imx662_table_fig('Table energy resolution multi', 'Table 4-4 equivalent: Energy resolution for multi-pixel events (n1234, 1-4 pixels)', H, D4, capR);

    % ---------- Tables 4-5 / 4-6: total counts; 4-7 / 4-8: count rates ----------
    DC5 = {}; DR = {};
    for i = 1:4
        L = getL(lines{i,1}, lines{i,2});
        DC5(i,:) = {nm(L.name), L.n1.area, ref.cnt1(i), L.n1234.area, ref.cnt4(i)}; %#ok<AGROW>
        DR(i,:)  = {nm(L.name), L.rate1, ref.rate1(i), L.rate1234, ref.rate4(i)}; %#ok<AGROW>
    end
    Hc = {'Emission line', 'Total counts, single-pixel', 'Lupo (2025) single', 'Total counts, multi-pixel', 'Lupo (2025) multi'};
    writecell([Hc; DC5], fullfile(out, 'Table_05_06_total_counts.csv'));
    imx662_table_fig('Table total counts', 'Tables 4-5/4-6 equivalent: Total counts under each emission line', Hc, DC5, ...
        ['Number of X-rays in each emission line = area under its Gaussian fit (Eq. 3.2: A = amplitude x sigma x sqrt(2 pi)). Totals depend on how long each source was measured ' ...
         sprintf('(this data: Fe-55 %.1f h, Po-210 %.1f h), so compare ratios between lines rather than absolute numbers with the thesis.', Fe.tTot/3600, Po.tTot/3600)]);
    Hr = {'Emission line', 'Count rate, single [counts/s/mm^2]', 'Lupo (2025) single', 'Count rate, multi [counts/s/mm^2]', 'Lupo (2025) multi'};
    writecell([Hr; DR], fullfile(out, 'Table_07_08_count_rates.csv'));
    imx662_table_fig('Table count rates', 'Tables 4-7/4-8 equivalent: Average count rate per emission line', Hr, DR, ...
        ['Average number of X-rays per second per square millimetre of sensor in each line (events within 3 sigma of the line, hot pixels removed, optical area 17.44 mm^2). ' ...
         'This is independent of measurement length, so it can be compared directly with the thesis and between sources.']);

    % ---------- Table 4-9: readout noise ----------
    DN = {'Average', RN.avg, 3.66; 'RMS', RN.rms, 3.77; 'Mode', RN.mode, 3.54; 'Median', RN.median, NaN};
    writecell([{'Parameter', 'Readout noise [e- RMS]', 'Lupo (2025) [e- RMS]'}; DN], fullfile(out, 'Table_09_readout_noise.csv'));
    imx662_table_fig('Table readout noise', 'Table 4-9 equivalent: Readout noise summary', {'Parameter', 'IMX662 readout noise [e- RMS]', 'Lupo (2025) [e- RMS]'}, DN, ...
        ['Readout noise of the pixels in electrons. Average = mean over pixels; RMS = root-mean-square (more sensitive to the few noisy pixels); mode = most common value (the typical pixel, ' ...
         'used in the energy-resolution budget). Converted from ADU with the Fe-55 gain and 3.65 eV per electron.']);

    % ---------- Gain, linearity, dark current summary ----------
    GV = Fe.GV;
    DS = {'Conversion gain from Mn K-alpha [ADU/eV]', Fe.gain, '0.2655'; ...
          'Conversion gain from X-ray linearity slope [ADU/eV]', XL.fit(1), '0.26482'; ...
          'Conversion gain from optical photon transfer [ADU/eV] (extra)', OL.gainPTC, '-'; ...
          'Pixel-to-pixel gain variation, Mn K-alpha [%]', GV.pct, '0.32'; ...
          'Gain-variation broadening sigma at 5900 eV [eV]', GV.sigmaEV, '18.88'; ...
          'Pixels used for gain variation', GV.nPix, '3600'; ...
          'X-ray linearity error [%]', XL.linErr, '0.301'; ...
          'Optical linearity error, 5-95% saturation [%]', OL.E1, '14.46'; ...
          'Optical linearity error, < 2000 ADU [%]', OL.E2, '0.22'; ...
          'Optical linearity error, 2000-4000 ADU [%]', OL.E3, '0.46'; ...
          'Optical saturation level [ADU]', OL.sat, '~4000'; ...
          'Dark current, lab method, full detector [e-/pixel/s]', DC.fits(1,1), '0.03 (not measurable)'; ...
          'Dark current, lab method, center region [e-/pixel/s]', DC.fits(2,1), 'not measurable'; ...
          'Dark current, lab method, edge columns [e-/pixel/s]', DC.fits(3,1), 'not measurable'; ...
          'Dark current, shot-noise method [e-/pixel/s] (extra)', DC.DCshot, '-'; ...
          'Hot dark-current pixels (> threshold) (extra)', DC.nHot, '-'};
    writecell([{'Quantity', 'Carissma (2026)', 'Lupo (2025) IMX662'}; DS], fullfile(out, 'Table_summary_gain_linearity_dark.csv'));
    imx662_table_fig('Table summary', 'Summary: gain, gain variation, linearity and dark current', {'Quantity', 'Carissma (2026)', 'Lupo (2025) IMX662'}, DS, ...
        ['Conversion gain: how many ADU (digital units) one eV of X-ray energy produces. Gain variation: how much this differs pixel to pixel. Linearity error (Eq. 3.4): how far the response ' ...
         'departs from a straight line, as a % of the range. Dark current: thermal signal per pixel per second with no light, at 25 C.']);

    % ---------- Tables 4-11, 4-13, 4-14, 4-15: resolution budget ----------
    B = XL.budget; DB = {};
    for i = 1:height(B)
        k = find(strcmp({'C K-alpha', 'F K-alpha', 'Mn K-alpha', 'Mn K-beta'}, B.Line{i}), 1);
        if strcmp(B.Events{i}, 'n1'), r = ref.un1; else, r = ref.un4; end
        rv = NaN; if ~isempty(k), rv = r(k); end
        DB(i,:) = {B.Line{i}, B.Events{i}, B.Measured_FWHM_eV(i), B.Fano_FWHM_eV(i), B.Calc_FWHM_Fano_Read_eV(i), B.Calc_FWHM_Fano_Read_Gain_eV(i), B.Unaccounted_sigma_eV(i), rv}; %#ok<AGROW>
    end
    Hb = {'Emission line', 'Events', 'Measured FWHM [eV]', 'Fano limit FWHM [eV]', sprintf('Calc. Fano+read FWHM [eV]'), 'Calc. Fano+read+gain [eV]', 'Unaccounted sigma [eV]', 'Lupo (2025) unacc. [eV]'};
    writecell([Hb; DB], fullfile(out, 'Table_10_15_energy_resolution_budget.csv'));
    imx662_table_fig('Table resolution budget', 'Tables 4-11/4-13/4-14/4-15 equivalent: What limits the energy resolution', Hb, DB, ...
        sprintf(['Calculated FWHM = 2.3548 x sqrt(sigma_Fano^2 + N x sigma_read^2 [+ sigma_gain^2]) (Eq. 3.5), with sigma_read = readout-noise mode (%.2f e- = %.1f eV) and N = mean pixels per event. ' ...
         '"Unaccounted sigma" = extra noise that must be added in quadrature to reach the measured width = sqrt(measured^2 - calculated^2)/2.3548.'], RN.mode, XL.sigReadEV));

    % ---------- Dark signal vs exposure ----------
    DD = num2cell([DC.texp DC.sig(:,1) DC.spat(:,1) DC.temp DC.sig(:,2) DC.sig(:,3)]);
    Hd = {'Exposure [s]', 'Dark signal, full [e-/pixel]', 'Spatial error [e-]', 'Temporal error [e-]', 'Dark signal, center [e-/pixel]', 'Dark signal, edges [e-/pixel]'};
    writecell([Hd; DD], fullfile(out, 'Table_dark_signal_vs_exposure.csv'));
    imx662_table_fig('Table dark signal', 'Dark-current measurement values (data behind the dark-current plots)', Hd, DD, ...
        ['Bias-subtracted, outlier-clipped mean dark signal per pixel at each exposure (the 15 us row is the bias reference itself, so it is zero). Spatial error = pixel-to-pixel RMS; ' ...
         'temporal error = frame-to-frame standard error (std / sqrt(25)). Center = outer 375-pixel columns removed; edges = only those columns.']);

    imx662_save_figs('7_Tables', out);
end
