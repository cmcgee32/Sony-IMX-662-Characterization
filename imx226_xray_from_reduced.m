function R = imx226_xray_from_reduced(folder, src, gainIn)
%IMX226_XRAY_FROM_REDUCED  IMX226 X-ray response (Fe-55 or Po-210) from reduced event data.
%
%   R = imx226_xray_from_reduced(folder, 'Po210')          % gain from F K-alpha (lab method for IMX226)
%   R = imx226_xray_from_reduced(folder, 'Fe55', gain)     % uses the Po-210 gain
%
%   The IMX226 event lists on Box are ~20 GB (Po-210) and ~50 GB (Fe-55), too
%   large to download, so they were streamed straight from Box and reduced in
%   one pass (imx226_box_reducer.js). Every event of every grades.fits file
%   used by the lab was read: Po-210 runs 2-353, Fe-55 data cubes 2-41 (see
%   PsyduckDataNotes.xlsx). Hot pixels are removed with the lab's own
%   per-event hot-pixel flags (hot_pixels_*.fits in the same zips), so the
%   cut is identical to the one used by Lupo (2025).
%
%   The reduced folder holds, as little-endian binary arrays:
%     S1/S4/S9   spectra of n1 / n1234 / all (<=9 pixel) events, 1-ADU bins, ALLSPPH
%     SH, SH1    spectra of hot-pixel events (all, n1)
%     np         event-size histogram (number of pixels above the split threshold)
%     lc         events per frame inside each line window, [nWindows x 1500 x 400 runs]
%     mapK       events per 4x4-pixel cluster inside window K (1000 x 750)
%     hot4       hot pixels per 4x4 cluster;  hc / hcGood  histograms of events per pixel
%     gv_*       Fe-55 only: mean single-pixel Mn K-alpha PH of every pixel with > 100 events
%   Output R has the same fields as imx662_xray_events, so imx662_plot_xray(R)
%   draws the lab-standard figures.
    c = imx662_cfg();
    isFe = strcmpi(src, 'Fe55');
    M = jsondecode(fileread(fullfile(folder, 'meta.json')));
    rd = @(n, t) fread_all(fullfile(folder, n), t);

    cfg.nx = 4000; cfg.ny = 3000; cfg.pixMM = 1.85e-3; cfg.areaMM2 = 7.4*5.55;   % Lupo (2025) Table 2-1
    cfg.frameTime = 0.6; cfg.cluster = 4;   % EXPTIME = 0.6 s in the FITS headers (time per frame used by the lab); 4x4 clusters
    nx = cfg.nx; ny = cfg.ny;

    maxADU = 4000;
    S1 = rd('S1.f64', 'double'); S4 = rd('S4.f64', 'double'); S9 = rd('S9.f64', 'double'); SH = rd('SH.f64', 'double');
    S.n1 = S1(2:maxADU+1); S.n1234 = S4(2:maxADU+1); S.all = S9(2:maxADU+1); S.hot = SH(2:maxADU+1); S.big = sum(S9(maxADU+2:end));
    np = rd('np.f64', 'double'); npixHist = np(2:10);

    if isFe
        L = struct('name', {'Si K\alpha','Si escape (Mn K\alpha)','Mn K\alpha','Mn K\beta'}, 'short', {'SiKa','SiEsc','MnKa','MnKb'}, ...
                   'E', {c.E.SiKa, c.E.SiEsc, c.E.MnKa, c.E.MnKb}, 'lc', {false,false,true,true}, 'win', {[],[],1,2});
        xMaxEV = 7800;
    else
        L = struct('name', {'C K\alpha','O K\alpha','F K\alpha'}, 'short', {'C','O','F'}, ...
                   'E', {c.E.C, c.E.O, c.E.F}, 'lc', {true,false,true}, 'win', {1,2,3});
        xMaxEV = 1000;
    end
    adu = (1:maxADU)';
    if nargin < 3 || isempty(gainIn), gainIn = 0.3654; end
    gain = gainIn;
    for pass = 1:2
        for j = 1:numel(L)
            L(j).n1    = fitLine(adu, S.n1,    L(j).E * gain);
            L(j).n1234 = fitLine(adu, S.n1234, L(j).E * gain);
        end
        if ~isFe   % lab method for the IMX226: conversion gain from the F K-alpha line
            jF = find(strcmp({L.short}, 'F')); gain = L(jF).n1.mu / c.E.F;
        else
            break
        end
    end
    gainMn = NaN; if isFe, jK = find(strcmp({L.short}, 'MnKa')); gainMn = L(jK).n1.mu / c.E.MnKa; end

    % frames and exposure
    frames = double(rd('frames.u32', 'uint32')); runs = M.runsDone(:)';
    nFr = frames(runs + 1); tTot = sum(nFr) * cfg.frameTime;
    lc = reshape(double(rd('lc.u32', 'uint32')), M.lcShape(3), 1500, 400);

    hot4 = reshape(double(rd('hot4.u16', 'uint16')), 1000, 750)';
    fwhmFano = @(E) 2.3548*sqrt(c.fano*c.epsSi*E);
    eV = 1/gain;
    for j = 1:numel(L)
        for cls = {'n1','n1234'}
            q = L(j).(cls{1});
            q.centerEV = q.mu*eV; q.centerEVerr = q.muErr*eV; q.sigmaEV = q.sig*eV; q.sigmaEVerr = q.sigErr*eV;
            q.fwhmEV = 2.3548*q.sigmaEV; q.fwhmEVerr = 2.3548*q.sigmaEVerr;
            L(j).(cls{1}) = q;
        end
        L(j).fano = fwhmFano(L(j).E);
        L(j).meanNpix = sum((1:4)' .* npixHist(1:4)) / sum(npixHist(1:4));   % mean pixels per n1234 event
        L(j).rate1 = L(j).n1.area / tTot / cfg.areaMM2;          % Gaussian-area counts / time / area
        L(j).rate1234 = L(j).n1234.area / tTot / cfg.areaMM2;
        if ~isempty(L(j).win)
            w = L(j).win; series = [];
            for r = runs, series = [series; squeeze(lc(w, 1:frames(r+1), r+1))']; end %#ok<AGROW>
            L(j).lcAll = series; L(j).rateAll = sum(series) / tTot / cfg.areaMM2;
            L(j).mapC = reshape(double(rd(sprintf('map%d.u32', w), 'uint32')), 1000, 750)';
            L(j).winADU = M.win(w, :);
        end
    end

    % events-per-pixel histogram and hot pixels
    hcAll = double(rd('hc.u32', 'uint32'));
    if isFe
        xmaxH = find(cumsum(hcAll) >= 0.9999*sum(hcAll), 1); xc = 0:xmaxH; hc = hcAll(1:xmaxH+1)';
        [~, im] = max(hc .* (xc > 0));
        [pg, pge] = imx662_gfit(xc(xc > 0), hc(xc > 0), [max(hc) xc(im) xc(im)/3]);
        hotThr = pg(2) + 3*pg(3); nHotGauss = sum(hcAll(xc(end)+1:end)) + sum(hc(xc > hotThr));
    else
        xc = 0:200; hc = hcAll(1:201)'; pg = [NaN NaN NaN]; pge = pg; hotThr = NaN; nHotGauss = NaN;
    end
    nHot = M.nHotPix;     % pixels with at least one lab-flagged hot event

    GV = struct();
    if isFe && exist(fullfile(folder, 'gv_meanph.f32'), 'file')
        mph = double(rd('gv_meanph.f32', 'single'));
        pixGain = mph / c.E.MnKa;
        be = linspace(prctile0(pixGain, 0.1), prctile0(pixGain, 99.9), 120);
        hg = histcounts(pixGain, be); bc = (be(1:end-1)+be(2:end))/2; [~, im] = max(hg);
        [pgv, pgve] = imx662_gfit(bc, hg, [max(hg) bc(im) std(pixGain)]);
        GV = struct('nPix', numel(pixGain), 'pixGain', pixGain, 'center', pgv(2), 'centerErr', pgve(2), 'sigma', pgv(3), ...
                    'pct', 100*pgv(3)/pgv(2), 'sigmaEV', pgv(3)/pgv(2)*c.E.MnKa, 'be', be, 'hg', hg, 'fit', pgv);
    end

    % IMX226 Fe-55: too few events per pixel for the per-pixel (> 100 events)
    % method, so the gain variation is estimated statistically from the spread
    % of per-pixel mean Mn K-alpha pulse heights minus counting noise.
    if isFe && ~isfield(GV, 'pct') && exist(fullfile(folder, 'gv_stat.json'), 'file')
        G = jsondecode(fileread(fullfile(folder, 'gv_stat.json'))); fn = fieldnames(G.res); r = G.res.(fn{1});
        GV = struct('pct', r.gainVarPct, 'sigmaEV', r.sigmaEV, 'nPix', r.pixels, ...
                    'method', 'statistical: spread of per-pixel mean Mn K-alpha PH minus counting noise, pixels with >= 10 events');
        fprintf('IMX226 Fe-55 gain variation (statistical, %d pixels): %.3f %% (%.1f eV at 5900 eV)\n', r.pixels, r.gainVarPct, r.sigmaEV);
    end

    R = struct('sensor', 'IMX226', 'cfg', cfg, 'src', src, 'gain', gain, 'gainMn', gainMn, 'L', L, 'S', S, ...
               'hotThr', hotThr, 'nHot', nHot, 'nHotGauss', nHotGauss, 'hotPct', 100*nHot/(nx*ny), 'hotC', hot4, ...
               'pg', pg, 'pge', pge, 'nFrames', sum(nFr), 'tTot', tTot, 'present', true(sum(nFr), 1), ...
               'npixHist', npixHist, 'GV', GV, 'xc', xc, 'hc', hc, 'xMaxEV', xMaxEV, 'isFe', isFe, ...
               'nEvents', M.nev, 'nHotEvents', M.nhot, 'runs', runs);
    fprintf('IMX226 %s: %d runs, %.2f h, %.3g events (%.1f%% flagged hot), gain %.5f ADU/eV, hot pixels %d (%.2f%%)\n', ...
            src, numel(runs), tTot/3600, M.nev, 100*M.nhot/M.nev, gain, nHot, R.hotPct);
end

function v = fread_all(fn, t)
    fid = fopen(fn, 'r', 'ieee-le'); v = fread(fid, Inf, ['*' t]); fclose(fid);
end

function q = fitLine(adu, h, mu0, ~)
    sig0 = 4 + 0.008*mu0;
    w = abs(adu - mu0) <= max(0.06*mu0, 3*sig0);
    hh = h; hh(~w) = 0; [pk, im] = max(hh); mu = adu(im); sig = sig0;
    if ~w(max(im-1,1)) || ~w(min(im+1,numel(w))), mu = mu0; pk = h(round(mu0)); end   % maximum at the window edge: no peak, start at the expected energy
    p = [pk mu sig 0]; perr = nan(1,4);
    for it = 1:4
        win = adu >= mu - 2*sig & adu <= mu + 2*sig;
        if nnz(win) < 6, win = adu >= mu - 6 & adu <= mu + 6; end
        try
            [pn, pen] = imx662_gfit(adu(win), h(win), [p(1) mu sig p(4)], true);
        catch
            break
        end
        % keep the fit on its own line: a weak line (O K-alpha) sitting on the
        % tail of a strong one (F K-alpha) must not drift onto the strong peak
        if abs(pn(2) - mu0) > max(0.06*mu0, 3*sig0) || abs(pn(3)) > 3*sig0, break; end
        p = pn; perr = pen;
        mu = p(2); sig = max(abs(p(3)), 1.5);
    end
    % a line that cannot be separated (no peak of its own near the expected
    % energy, e.g. O K-alpha under the F K-alpha tail in n1234) is reported as NaN
    if abs(p(2) - mu0) > max(0.04*mu0, 2*sig0) || any(~isfinite(perr(2:3)))
        p(:) = NaN; perr(:) = NaN;
    end
    q.A = p(1); q.mu = p(2); q.sig = abs(p(3)); q.bkg = p(4);
    q.muErr = perr(2); q.sigErr = perr(3);
    q.area = q.A * q.sig * sqrt(2*pi);
    q.areaErr = q.area * sqrt((perr(1)/max(q.A,eps))^2 + (perr(3)/max(q.sig,eps))^2);
    q.e33 = q.sig / sqrt(max(q.area, 1));
end

function v = prctile0(x, p)
    x = sort(x(:)); v = x(max(1, min(numel(x), round(p/100*numel(x)))));
end
