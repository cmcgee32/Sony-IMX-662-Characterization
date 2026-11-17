function R = imx662_xray_events(src, gainIn)
%IMX662_XRAY_EVENTS  X-ray response of the IMX662 from event-mode data (Fe-55 or Po-210).
%
%   R = imx662_xray_events('Fe55')            % fits Mn K-alpha -> conversion gain
%   R = imx662_xray_events('Po210', gain)     % uses the Fe-55 gain (lab method)
%
% Follows the lab method (Lupo 2025, Sec. 3.1-3.3):
%   * Hot pixels: Fe-55 -> histogram of events per pixel, Gaussian fit,
%     threshold = center + 3 sigma.  Po-210 -> any pixel with >= 10 events.
%   * Event classes by number of pixels above the 9 ADU split threshold:
%     n1 = single-pixel, n1234 = 1-4 pixels, all = 1-9 pixels.
%     Pulse height = all-split pulse height (ALLSPPH) [ADU].
%   * Gaussian fits to each emission line -> center, sigma, FWHM, area.
%   * Conversion gain [ADU/eV] = Mn K-alpha center [ADU] / 5900 eV.
%   * Light curves: events per frame within +-3 sigma of each line,
%     / 0.6 s / optical area -> counts/s/mm^2, rebinned to 0.1% of the run.
%   * Count-rate heat maps averaged over 3x3-pixel clusters.
%   * Gain variance (Fe-55, Mn K-alpha single-pixel events): per-pixel gain
%     for every pixel with > 100 events.
%
% Reads the zip directly; one file at a time is extracted to a temp folder.

    c = imx662_cfg();
    isFe = strcmpi(src, 'Fe55');
    if isFe, zipPath = c.zipFe; else, zipPath = c.zipPo; end
    if nargin < 2 || isempty(gainIn), gainIn = c.gainRef; end
    tmp = fullfile(tempdir, ['imx662_ev_' src]); if exist(tmp,'dir'), rmdir(tmp,'s'); end; mkdir(tmp);
    cln = onCleanup(@() rmdir(tmp, 's'));
    nx = c.nx; ny = c.ny; Npix = nx*ny; T = c.frameTime;

    % Emission lines to fit
    if isFe
        L = struct('name', {'Si K\alpha','Si escape (Mn K\alpha)','Mn K\alpha','Mn K\beta'}, ...
                   'short', {'SiKa','SiEsc','MnKa','MnKb'}, ...
                   'E', {c.E.SiKa, c.E.SiEsc, c.E.MnKa, c.E.MnKb}, 'lc', {false,false,true,true});
        xMaxEV = 7800;
    else
        L = struct('name', {'C K\alpha','O K\alpha','F K\alpha'}, 'short', {'C','O','F'}, ...
                   'E', {c.E.C, c.E.O, c.E.F}, 'lc', {true,false,true});
        xMaxEV = 1000;
    end

    names = imx662_zip_list(zipPath, '_grades\.fits$');
    names = names(~contains(names, 'OlderFiles'));
    idx = cellfun(@(s) str2double(regexp(s, '_(\d+)_grades\.fits$', 'tokens', 'once')), names);
    [idx, o] = sort(idx); names = names(o); nC = numel(names);
    fprintf('\n=== %s: %d chunks ===\n', src, nC);

    % ---------------- PASS 1: events per pixel (hot pixels), frame range ----------------
    cnt = zeros(Npix, 1); fmin = inf; fmax = -inf; chunkFr = zeros(nC, 2); t0 = tic;
    for k = 1:nC
        Tb = readCols(zipPath, names{k}, tmp, [1 2 3]);
        fr = double(Tb{1}); li = double(Tb{2})*ny + double(Tb{3}) + 1;
        li = li(li >= 1 & li <= Npix);
        cnt = cnt + accumarray(li, 1, [Npix 1]);
        if ~isempty(fr), chunkFr(k,:) = [min(fr) max(fr)]; fmin = min(fmin, min(fr)); fmax = max(fmax, max(fr)); end
        if mod(k, 10) == 1 || k == nC, fprintf('  pass 1/3  chunk %d/%d  (%.0f s)\n', k, nC, toc(t0)); end
    end
    % Compact frame index: chunks are laid end to end (frame counters can jump between runs)
    span = max(min(chunkFr(:,2) - chunkFr(:,1) + 1, 3000), 0);
    span(~any(chunkFr, 2) & span <= 1) = 0;
    off = [0; cumsum(span(1:end-1))];
    nF = sum(span); present = true(nF, 1);
    nFrames = nF; tTot = nFrames * T;
    fprintf('Frames: %d (%.2f h of exposure); raw frame counter range %d ... %d\n', nF, tTot/3600, fmin, fmax);

    % Hot-pixel threshold
    if isFe
        e = 0:1:max(cnt); hc = histcounts(cnt, [e e(end)+1]); xc = e;
        [~, im] = max(hc .* (xc > 0));
        [pg, pge] = imx662_gfit(xc(xc > 0), hc(xc > 0), [max(hc) im std(cnt(cnt>0))]);
        hotThr = pg(2) + 3*pg(3);
        hot = cnt > hotThr;
    else
        hotThr = 10; hot = cnt >= hotThr; pg = [NaN NaN NaN]; pge = pg; xc = 0:max(min(max(cnt), 200), 20); hc = histcounts(cnt, [xc xc(end)+1]);
    end
    nHot = nnz(hot);
    fprintf('Hot pixels: threshold %.1f events/pixel -> %d pixels (%.2f%%)\n', hotThr, nHot, 100*nHot/Npix);

    % ---------------- PASS 2: spectra ----------------
    maxADU = 4000; if ~isFe, maxADU = 4000; end
    S.n1 = zeros(maxADU,1); S.n1234 = S.n1; S.all = S.n1; S.hot = S.n1; S.big = 0;
    npixHist = zeros(9,1); t0 = tic;
    for k = 1:nC
        Tb = readCols(zipPath, names{k}, tmp, [2 3 5 8]);
        li = double(Tb{1})*ny + double(Tb{2}) + 1; np = double(Tb{3}); ph = round(double(Tb{4}));
        ok = li >= 1 & li <= Npix & ph >= 1;
        li = li(ok); np = np(ok); ph = ph(ok);
        isHot = hot(li);
        big = ph > maxADU; S.big = S.big + nnz(big & ~isHot); ph = min(ph, maxADU);
        g = ~isHot;
        S.n1    = S.n1    + accumarray(ph(g & np == 1), 1, [maxADU 1]);
        S.n1234 = S.n1234 + accumarray(ph(g & np <= 4), 1, [maxADU 1]);
        S.all   = S.all   + accumarray(ph(g & np <= 9), 1, [maxADU 1]);
        S.hot   = S.hot   + accumarray(ph(isHot & np <= 9), 1, [maxADU 1]);
        npixHist = npixHist + accumarray(min(max(np(g),1),9), 1, [9 1]);
        if mod(k, 10) == 1 || k == nC, fprintf('  pass 2/3  chunk %d/%d  (%.0f s)\n', k, nC, toc(t0)); end
    end
    adu = (1:maxADU)';

    % ---------------- Line fits (in ADU), then gain ----------------
    gain = gainIn;
    for pass = 1:2
        for j = 1:numel(L)
            L(j).n1    = fitLine(adu, S.n1,    L(j).E * gain);
            L(j).n1234 = fitLine(adu, S.n1234, L(j).E * gain);
        end
        if isFe
            jK = find(strcmp({L.short}, 'MnKa'));
            gain = L(jK).n1.mu / c.E.MnKa;            % lab definition: center / known energy
        else
            break
        end
    end
    fprintf('Conversion gain used: %.5f ADU/eV  (lab reference %.4f)\n', gain, c.gainRef);
    eVperADU = 1/gain;

    % ---------------- PASS 3: light curves, heat maps, per-pixel gain ----------------
    nL = numel(L);
    for j = 1:nL
        w = L(j).n1234; L(j).win = [w.mu - 3*w.sig, w.mu + 3*w.sig];
        L(j).lcAll = zeros(nF,1); L(j).lc1 = zeros(nF,1); L(j).lc1234 = zeros(nF,1);
        L(j).map = zeros(Npix,1); L(j).npixSum = 0; L(j).npixN = 0;
    end
    gSum = zeros(Npix,1); gN = zeros(Npix,1); jK = find(strcmp({L.short}, 'MnKa'));
    t0 = tic;
    for k = 1:nC
        Tb = readCols(zipPath, names{k}, tmp, [1 2 3 5 8]);
        fi = double(Tb{1}) - chunkFr(k,1) + 1; fi(fi > span(k)) = 0; fi = fi + off(k) .* (fi > 0);
        li = double(Tb{2})*ny + double(Tb{3}) + 1;
        np = double(Tb{4}); ph = double(Tb{5});
        ok = li >= 1 & li <= Npix & fi >= 1 & fi <= nF;
        fi = fi(ok); li = li(ok); np = np(ok); ph = ph(ok);
        g = ~hot(li) & np <= 9;
        for j = 1:nL
            inW = g & ph >= L(j).win(1) & ph <= L(j).win(2);
            L(j).lcAll  = L(j).lcAll  + accumarray(fi(inW), 1, [nF 1]);
            L(j).lc1    = L(j).lc1    + accumarray(fi(inW & np == 1), 1, [nF 1]);
            L(j).lc1234 = L(j).lc1234 + accumarray(fi(inW & np <= 4), 1, [nF 1]);
            L(j).map    = L(j).map    + accumarray(li(inW), 1, [Npix 1]);
            L(j).npixSum = L(j).npixSum + sum(np(inW & np <= 4)); L(j).npixN = L(j).npixN + nnz(inW & np <= 4);
        end
        if isFe
            w1 = L(jK).n1; sel = g & np == 1 & abs(ph - w1.mu) <= 3*w1.sig;
            gSum = gSum + accumarray(li(sel), ph(sel), [Npix 1]);
            gN   = gN   + accumarray(li(sel), 1, [Npix 1]);
        end
        if mod(k, 10) == 1 || k == nC, fprintf('  pass 3/3  chunk %d/%d  (%.0f s)\n', k, nC, toc(t0)); end
    end

    % ---------------- Tables ----------------
    fwhmFano = @(E) 2.3548*sqrt(c.fano*c.epsSi*E);
    for j = 1:nL
        for cls = {'n1','n1234'}
            q = L(j).(cls{1});
            q.centerEV = q.mu * eVperADU;   q.centerEVerr = q.muErr * eVperADU;
            q.sigmaEV  = q.sig * eVperADU;  q.sigmaEVerr  = q.sigErr * eVperADU;
            q.fwhmEV   = 2.3548*q.sigmaEV;  q.fwhmEVerr   = 2.3548*q.sigmaEVerr;
            L(j).(cls{1}) = q;
        end
        L(j).fano = fwhmFano(L(j).E);
        L(j).rate1    = sum(L(j).lc1(present))    / tTot / c.areaMM2;
        L(j).rate1234 = sum(L(j).lc1234(present)) / tTot / c.areaMM2;
        L(j).rateAll  = sum(L(j).lcAll(present))  / tTot / c.areaMM2;
        L(j).meanNpix = L(j).npixSum / max(L(j).npixN, 1);
    end
    fprintf('\n%-24s %16s %16s %10s %14s %14s %12s %12s\n', 'Line', 'FWHM n1 [eV]', 'FWHM n1234 [eV]', 'Fano [eV]', ...
            'Counts n1', 'Counts n1234', 'Rate n1', 'Rate n1234');
    for j = 1:nL
        fprintf('%-24s %8.1f +- %5.1f %8.1f +- %5.1f %10.2f %14.0f %14.0f %12.4f %12.4f\n', strrep(L(j).name,'\',''), ...
            L(j).n1.fwhmEV, L(j).n1.fwhmEVerr, L(j).n1234.fwhmEV, L(j).n1234.fwhmEVerr, L(j).fano, ...
            L(j).n1.area, L(j).n1234.area, L(j).rate1, L(j).rate1234);
    end

    % Gain variance (Fe-55)
    GV = struct();
    if isFe
        q = gN > 100;
        pixGain = (gSum(q) ./ gN(q)) / c.E.MnKa;      % ADU/eV per pixel
        GV.nPix = nnz(q); GV.pixGain = pixGain;
        be = linspace(prctile0(pixGain, 0.1), prctile0(pixGain, 99.9), 120);
        hg = histcounts(pixGain, be); bc = (be(1:end-1)+be(2:end))/2;
        [~, im] = max(hg);
        [pgv, pgve] = imx662_gfit(bc, hg, [max(hg) bc(im) std(pixGain)]);
        GV.center = pgv(2); GV.centerErr = pgve(2); GV.sigma = pgv(3);
        GV.pct = 100 * GV.sigma / GV.center;
        GV.sigmaEV = GV.pct/100 * c.E.MnKa;
        GV.be = be; GV.hg = hg; GV.fit = pgv;
        fprintf('Gain variance: %d pixels with >100 Mn K-alpha single-pixel events; center %.5f ADU/eV, sigma %.5f -> %.3f%% -> %.1f eV\n', ...
                GV.nPix, GV.center, GV.sigma, GV.pct, GV.sigmaEV);
    end

    R = struct('src', src, 'gain', gain, 'L', L, 'S', S, 'hotThr', hotThr, 'nHot', nHot, 'hotPct', 100*nHot/Npix, ...
               'cnt', reshape(cnt, ny, nx), 'hotMap', reshape(hot, ny, nx), 'pg', pg, 'pge', pge, ...
               'nFrames', nFrames, 'tTot', tTot, 'present', present, 'npixHist', npixHist, 'GV', GV, ...
               'xc', xc, 'hc', hc, 'xMaxEV', xMaxEV, 'isFe', isFe);
    % Figures are made separately: imx662_plot_xray(R)
end

% ======================================================================
function q = fitLine(adu, h, mu0, ~)
    sig0 = 4 + 0.008*mu0;
    w = abs(adu - mu0) <= max(0.06*mu0, 3*sig0);
    hh = h; hh(~w) = 0; [pk, im] = max(hh); mu = adu(im); sig = sig0;
    p = [pk mu sig 0]; perr = nan(1,4);
    for it = 1:4
        win = adu >= mu - 2*sig & adu <= mu + 2*sig;
        if nnz(win) < 6, win = adu >= mu - 6 & adu <= mu + 6; end
        try
            [p, perr] = imx662_gfit(adu(win), h(win), [p(1) mu sig p(4)], true);
        catch
            break
        end
        mu = p(2); sig = max(abs(p(3)), 1.5);
    end
    q.A = p(1); q.mu = p(2); q.sig = abs(p(3)); q.bkg = p(4);
    q.muErr = perr(2); q.sigErr = perr(3);
    q.area = q.A * q.sig * sqrt(2*pi);                   % Eq. 3.2 (bins are 1 ADU wide)
    q.areaErr = q.area * sqrt((perr(1)/max(q.A,eps))^2 + (perr(3)/max(q.sig,eps))^2);
    q.e33 = q.sig / sqrt(max(q.area, 1));                % Eq. 3.3 (weight for X-ray linearity)
end

function Tb = readCols(zipPath, name, tmp, cols)
    f = imx662_zip_extract(zipPath, name, tmp);
    cl = onCleanup(@() delete(f));
    Tb = fitsread(f, 'binarytable', 1, 'TableColumns', cols);
end

function v = prctile0(x, p)
    x = sort(x(:)); v = x(max(1, min(numel(x), round(p/100*numel(x)))));
end

