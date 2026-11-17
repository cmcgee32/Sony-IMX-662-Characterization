function R = imx662_xray_spectrum(zipPath, label, biasEvery)
%IMX662_XRAY_SPECTRUM  Event-mode spectrum analysis (Fe-55 or Po-210 data).
%
%   R = imx662_xray_spectrum(zipPath, 'Fe55')   or   (..., 'Po210')
%
% What the data are
%   The camera ran in event mode: every 0.6 s frame was bias-subtracted, and
%   each local maximum above EVTTH = 30 DN became an "event" with a 3x3 pixel
%   island around it. Per chunk of 1500 frames (~5 min) it saved
%     *_grades.fits : FRAME, X, Y, EVGEOM, NUMPIX, ASCA grade, PH, ALLSPPH, PHAS(3x3)
%     *_evts.fits   : same events, raw 3x3 only (not needed here)
%     *_biasmap.fits: per-pixel BIAS and NOISE maps used for that chunk
%   PH = summed charge of the pixels above the split threshold (9 DN) [DN].
%   ASCA grade: 0 = single pixel, 2/3/4 = 2-pixel splits, 6 = L/square
%   (all "good" X-ray grades); 1, 5, 7 = diagonal/large = background, cosmics,
%   or very extended charge clouds (e.g. alpha particles).
%
% What this computes
%   pulse-height spectra by grade, event rate vs time, where events land
%   on the sensor, cluster sizes, and bias/noise stability over the run.
%   For Fe-55 it fits the Mn K-alpha (5.895 keV) and K-beta (6.490 keV)
%   lines in the single-pixel spectrum, giving the conversion gain
%   K = (5895 eV / 3.65 eV per e-) / peak  [e-/DN]  and the energy resolution.

    if nargin < 3, biasEvery = 1; end
    tmp = fullfile(tempdir, ['imx662_' label]); if exist(tmp,'dir'), rmdir(tmp,'s'); end; mkdir(tmp);
    cleanupTmp = onCleanup(@() rmdir(tmp, 's'));

    gnames = imx662_zip_list(zipPath, '_grades\.fits$');
    gnames = gnames(~contains(gnames, 'OlderFiles'));
    idx = cellfun(@(s) str2double(regexp(s, '_(\d+)_grades\.fits$', 'tokens', 'once')), gnames);
    [idx, o] = sort(idx); gnames = gnames(o);
    nC = numel(gnames);
    fprintf('%s: %d chunks of event data\n', label, nC);

    maxPH = 40000; edges = 0:1:maxPH;                 % 1-DN bins
    Hg0 = zeros(1, maxPH); Hg234 = Hg0; Hg6 = Hg0; Hbad = Hg0;
    Hnpix = zeros(1, 10);
    cmap = zeros(1080, 1920);                         % where events land (good grades)
    cmapBad = cmap;
    nEv = zeros(nC,1); nFr = nEv; nGood = nEv; biasMed = nan(nC,1); noiseMed = biasMed;
    t0 = tic;
    for k = 1:nC
        f = imx662_zip_extract(zipPath, gnames{k}, tmp);
        try
            T = fitsread(f, 'binarytable', 1, 'TableColumns', [1 2 3 5 6 7]);
        catch ME
            delete(f); warning('%s: %s', gnames{k}, ME.message); continue;
        end
        delete(f);
        fr = double(T{1}); x = double(T{2}); y = double(T{3});
        npix = double(T{4}); g = double(T{5}); ph = double(T{6});
        nEv(k) = numel(ph);
        if ~isempty(fr), nFr(k) = max(fr) - min(fr) + 1; end
        good = ismember(g, [0 2 3 4 6]); nGood(k) = nnz(good);
        phb = min(max(round(ph), 1), maxPH);
        Hg0   = Hg0   + accumarray(phb(g==0), 1, [maxPH 1])';
        Hg234 = Hg234 + accumarray(phb(ismember(g,[2 3 4])), 1, [maxPH 1])';
        Hg6   = Hg6   + accumarray(phb(g==6), 1, [maxPH 1])';
        Hbad  = Hbad  + accumarray(phb(~good), 1, [maxPH 1])';
        Hnpix = Hnpix + accumarray(min(max(npix,1),10), 1, [10 1])';
        xi = min(max(x+1,1),1920); yi = min(max(y+1,1),1080);
        cmap    = cmap    + accumarray([yi(good) xi(good)], 1, [1080 1920]);
        cmapBad = cmapBad + accumarray([yi(~good) xi(~good)], 1, [1080 1920]);
        if mod(k-1, biasEvery) == 0
            bn = strrep(gnames{k}, '_grades.fits', '_biasmap.fits');
            try
                fb = imx662_zip_extract(zipPath, bn, tmp);
                B = fitsread(fb, 'image', 1); N = fitsread(fb, 'image', 2); delete(fb);
                biasMed(k) = median(B(:)); noiseMed(k) = median(N(:));
            catch
            end
        end
        fprintf('  chunk %3d/%d (#%d): %8d events, %5.1f%% good grades  (%.0f s elapsed)\n', ...
                k, nC, idx(k), nEv(k), 100*nGood(k)/max(nEv(k),1), toc(t0));
    end
    frameTime = 0.6;   % s (EXPTIME 600000 us)
    rate = nEv ./ (nFr * frameTime);
    phc = 1:maxPH;     % bin centres

    R = struct('label',label,'chunk',idx,'nEv',nEv,'nFr',nFr,'nGood',nGood,'rate',rate, ...
               'Hg0',Hg0,'Hg234',Hg234,'Hg6',Hg6,'Hbad',Hbad,'Hnpix',Hnpix, ...
               'cmap',cmap,'cmapBad',cmapBad,'biasMed',biasMed,'noiseMed',noiseMed);

    fprintf('\n===== %s RESULTS =====\n', label);
    fprintf('Total events: %d in %d frames (%.1f h); mean rate %.1f events/s\n', ...
            sum(nEv), sum(nFr), sum(nFr)*frameTime/3600, sum(nEv)/(sum(nFr)*frameTime));
    fprintf('Grade mix: g0 %.1f%%, g2-4 %.1f%%, g6 %.1f%%, other %.1f%%\n', ...
            100*[sum(Hg0) sum(Hg234) sum(Hg6) sum(Hbad)]/sum([Hg0 Hg234 Hg6 Hbad]));

    % ---------------- Fe-55 line fit ----------------
    if strcmpi(label, 'Fe55')
        [~, ipk] = max(Hg0 .* (phc > 300));          % K-alpha peak position
        w = phc > 0.85*ipk & phc < 1.25*ipk;
        xw = phc(w); yw = Hg0(w);
        % two Gaussians (K-alpha, K-beta at 6490/5895 x) + flat background
        model = @(p, x) p(1)*exp(-0.5*((x-p(2))/p(3)).^2) + ...
                        p(4)*exp(-0.5*((x-p(2)*p(6))/p(3)).^2) + p(5);
        p0 = [max(yw) ipk 0.03*ipk 0.15*max(yw) min(yw) 6490/5895];
        cost = @(p) sum((yw - model(p, xw)).^2 ./ max(yw, 1));
        p = fminsearch(cost, p0, optimset('MaxFunEvals', 2e4, 'MaxIter', 2e4));
        muA = p(2); sig = abs(p(3)); muB = p(2)*p(6);
        nA = 5895/3.65;                               % electrons per K-alpha photon
        K = nA / muA;                                 % e-/DN
        fwhm_eV = 2.3548 * sig * 5895 / muA;
        R.fit = p; R.K_eDN = K; R.KalphaDN = muA; R.KbetaDN = muB; R.FWHM_eV = fwhm_eV;
        fprintf('Mn K-alpha peak: %.1f DN   K-beta: %.1f DN   (ratio %.4f, expected %.4f)\n', ...
                muA, muB, muB/muA, 6490/5895);
        fprintf('Conversion gain from Fe-55: K = %.4f e-/DN  (1 DN = %.2f eV)\n', K, K*3.65);
        fprintf('Energy resolution (single-pixel) FWHM at 5.9 keV: %.0f eV (%.1f%%)\n', fwhm_eV, 100*fwhm_eV/5895);
        evPerDN = 5895 / muA;
    else
        evPerDN = NaN;
    end

    % ---------------- Figures ----------------
    figure('Name', [label ' spectrum'], 'Position', [40 40 1300 480]);
    if strcmpi(label, 'Fe55'), xmax = 2.0*R.KalphaDN; else, xmax = maxPH; end
    sm = @(h) movmean(h, 3);
    subplot(1,2,1);
    semilogy(phc, max(sm(Hg0),0.5), phc, max(sm(Hg234),0.5), phc, max(sm(Hg6),0.5), phc, max(sm(Hbad),0.5));
    grid on; xlim([0 xmax]);
    xlabel('Pulse height PH [DN]'); ylabel('Events per DN');
    legend('grade 0 (single pixel)','grades 2-4 (split)','grade 6','grades 1,5,7 (other)');
    title(sprintf('%s pulse-height spectrum by grade', label));
    subplot(1,2,2);
    if strcmpi(label, 'Fe55')
        plot(phc, Hg0, '.-'); hold on; grid on;
        plot(xw, model(p, xw), 'r-', 'LineWidth', 1.5);
        xlim([0.75 1.35]*R.KalphaDN);
        xlabel('PH [DN]'); ylabel('Events per DN');
        title(sprintf('Grade-0 Mn K\\alpha/K\\beta fit: K = %.3f e^-/DN, FWHM %.0f eV', R.K_eDN, R.FWHM_eV));
        ax2 = gca; ax2.XAxis.Exponent = 0;
    else
        loglog(phc, max(Hg0+Hg234+Hg6+Hbad, 0.5)); grid on;
        xlabel('Pulse height PH [DN] (log)'); ylabel('Events per DN');
        title('All grades, log-log: low-PH noise/background vs high-PH particle tail');
    end

    figure('Name', [label ' rate & stability'], 'Position', [70 70 1300 420]);
    tmin = cumsum(nFr) * frameTime / 60;
    subplot(1,3,1); plot(tmin, rate, '.-'); grid on;
    xlabel('Time [min]'); ylabel('Events / s'); title('Event rate per chunk');
    subplot(1,3,2); plot(tmin, 100*nGood./max(nEv,1), '.-'); grid on;
    xlabel('Time [min]'); ylabel('% good grades (0,2,3,4,6)'); title('Grade mix vs time');
    subplot(1,3,3); ok = isfinite(biasMed);
    yyaxis left;  plot(tmin(ok), biasMed(ok), 'o-'); ylabel('Median bias [DN]');
    yyaxis right; plot(tmin(ok), noiseMed(ok), 's-'); ylabel('Median noise [DN]');
    grid on; xlabel('Time [min]'); title('Bias & read-noise stability');

    figure('Name', [label ' spatial'], 'Position', [100 100 1300 420]);
    subplot(1,3,1); imagesc(imgaussfilt0(cmap, 8)); axis image off; colorbar;
    title('Good-grade events per pixel (smoothed)');
    subplot(1,3,2); imagesc(imgaussfilt0(cmapBad, 8)); axis image off; colorbar;
    title('Other-grade events per pixel (smoothed)');
    subplot(1,3,3); bar(1:10, Hnpix/sum(Hnpix)*100); grid on;
    xlabel('Pixels above split threshold (NUMPIX)'); ylabel('% of events');
    title('Event size distribution');
    colormap(parula);
    R.evPerDN = evPerDN;
end

function B = imgaussfilt0(A, s)
% Gaussian smoothing with base MATLAB only (no Image Processing Toolbox).
    h = exp(-0.5*((-3*s:3*s)/s).^2); h = h / sum(h);
    B = conv2(h, h, A, 'same');
end
