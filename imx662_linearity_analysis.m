function R = imx662_linearity_analysis(zipPath)
%IMX662_LINEARITY_ANALYSIS  Linearity + photon-transfer analysis of the IMX662.
%
%   R = imx662_linearity_analysis('/path/IMX662_Linearity_20250624.zip')
%
% What the data are
%   101 pairs of 25-frame cubes, "ShutterOpen" (light) and "ShutterClosed"
%   (dark), at exposure times from 2.0 s down to 16 us (the number at the end
%   of each file name is the exposure in microseconds; the header EXPTIME only
%   holds the session's starting value). Gain 9, same as all other datasets.
%
% What this computes, for every exposure
%   signal  S  = mean(open) - mean(closed)      in a central ROI      [DN]
%   shot-noise variance  V = var_temporal(open) - var_temporal(closed) [DN^2]
%   (temporal variance from consecutive-frame differences, so fixed-pattern
%    noise and slow light drifts cancel)
% and from those
%   * linearity: S vs t, residuals from a straight line (percent)
%   * full well / saturation: where S stops growing
%   * photon transfer curve (PTC): V vs S.  For a shot-noise-limited sensor
%     V = S / K, so the slope gives the conversion gain K [e-/DN]
%   * PRNU: pixel-to-pixel response non-uniformity at mid signal

    if nargin < 1
        zipPath = '/Users/carissmamcgee/Downloads/IMX662_Linearity_20250624.zip';
    end
    roiHalf = 256;                         % ROI = central 512 x 512 pixels
    tmp = fullfile(tempdir, 'imx662_lin'); if exist(tmp,'dir'), rmdir(tmp,'s'); end; mkdir(tmp);
    cleanupTmp = onCleanup(@() rmdir(tmp, 's'));

    names = imx662_zip_list(zipPath, 'datacube_Shutter(Open|Closed)_[\d.]+\.fits$');
    tok = regexp(names, 'datacube_Shutter(Open|Closed)_([\d.]+)\.fits$', 'tokens', 'once');
    isOpen = cellfun(@(c) strcmp(c{1}, 'Open'), tok);
    tus    = cellfun(@(c) str2double(c{2}), tok);
    texps  = unique(tus);
    n = numel(texps);
    fprintf('Linearity: %d exposures (%g us ... %g s)\n', n, min(texps), max(texps)/1e6);

    S = nan(n,1); So = S; Sc = S; Vo = S; Vc = S; satFrac = S; prnu = S; ffMean = S;
    keepMap = []; keepT = NaN;
    tMid = 0.3e6;                                       % exposure for the PRNU/illumination map
    [~, iMid] = min(abs(texps - tMid));

    t0 = tic;
    for k = 1:n
        io = find(isOpen  & tus == texps(k), 1);
        ic = find(~isOpen & tus == texps(k), 1);
        if isempty(io) || isempty(ic), continue; end
        [mO, vO, fO] = cubeStats(zipPath, names{io}, tmp, roiHalf);
        [mC, vC, fC] = cubeStats(zipPath, names{ic}, tmp, roiHalf);
        d = mO.roi - mC.roi;
        So(k) = mean(mO.roi(:));  Sc(k) = mean(mC.roi(:));
        S(k)  = mean(d(:));
        Vo(k) = vO;  Vc(k) = vC;
        satFrac(k) = mean(mO.roi(:) > 4000);
        prnu(k) = std(d(:)) / mean(d(:));
        ffMean(k) = fO - fC;
        if k == iMid, keepMap = mO.full - mC.full; keepT = texps(k); end
        fprintf('  %3d/%d  t = %9.6f s   S = %7.1f DN   V = %7.1f DN^2   (%.0f s elapsed)\n', ...
                k, n, texps(k)/1e6, S(k), Vo(k)-Vc(k), toc(t0));
    end

    t = texps/1e6;  V = Vo - Vc;
    ok = isfinite(S);
    % Linear range for fits: 2% .. 70% of the saturated level
    Smax = max(S(ok));
    lin  = ok & S > 0.02*Smax & S < 0.70*Smax & satFrac < 0.001;
    pl   = polyfit(t(lin), S(lin), 1);
    resid = 100 * (S - polyval(pl, t)) ./ polyval(pl, t);
    % PTC: fit V = S/K + c in the same range (shot-noise dominated part)
    pp   = polyfit(S(lin), V(lin), 1);
    K    = 1 / pp(1);                                   % e-/DN
    fullWellDN = Smax;
    readNoiseDN = sqrt(max(Vc(t == min(t)), 0)/1);      % from the shortest dark pair

    R = struct('t',t,'S',S,'So',So,'Sc',Sc,'V',V,'Vo',Vo,'Vc',Vc,'satFrac',satFrac, ...
               'prnu',prnu,'ffMean',ffMean,'linFit',pl,'resid',resid,'ptcFit',pp, ...
               'K_eDN',K,'fullWellDN',fullWellDN,'fullWell_e',fullWellDN*K, ...
               'readNoiseDN',readNoiseDN,'linMask',lin,'illumMap',keepMap,'illumT',keepT/1e6);

    fprintf('\n===== LINEARITY RESULTS =====\n');
    fprintf('Signal rate (linear fit): %.1f DN/s  (offset %.1f DN)\n', pl(1), pl(2));
    fprintf('Max nonlinearity in fit range: %.2f %%\n', max(abs(resid(lin))));
    fprintf('Saturation (max signal above dark): %.0f DN\n', Smax);
    fprintf('PTC conversion gain K = %.3f e-/DN  (=> full well ~ %.0f e-)\n', K, Smax*K);
    fprintf('PRNU at %.2f s: %.2f %%\n', keepT/1e6, 100*prnu(iMid));

    % ---------------- Figures ----------------
    figure('Name','Linearity','Position',[40 40 1300 450]);
    subplot(1,2,1);
    plot(t, S, 'o', 'MarkerSize', 4); hold on; grid on;
    plot(t(lin), S(lin), 'o', 'MarkerFaceColor', [0 .45 .74], 'MarkerSize', 4);
    tt = linspace(0, max(t), 200); plot(tt, polyval(pl, tt), 'r-');
    xlabel('Exposure [s]'); ylabel('Signal = open - closed [DN]');
    title(sprintf('Signal vs exposure: %.0f DN/s, saturates at ~%.0f DN', pl(1), Smax));
    legend('all','used for fit','linear fit','Location','southeast');
    subplot(1,2,2);
    plot(S(ok), resid(ok), 'o-', 'MarkerSize', 4); grid on; hold on;
    yline(0,'k-'); yline([-1 1],'k:');
    xlabel('Signal [DN]'); ylabel('Deviation from straight line [%]');
    title('Non-linearity (dotted = \pm1%)'); ylim([-10 10]);

    figure('Name','Photon transfer curve','Position',[70 70 1300 450]);
    subplot(1,2,1);
    loglog(S(ok & S>0), V(ok & S>0), 'o', 'MarkerSize', 4); hold on; grid on;
    ss = logspace(0, log10(Smax), 100); loglog(ss, ss/K, 'r-');
    xlabel('Signal [DN]'); ylabel('Shot-noise variance [DN^2]');
    title(sprintf('Photon transfer curve: K = %.2f e^-/DN', K));
    legend('data','V = S/K (slope 1)','Location','northwest');
    subplot(1,2,2);
    plot(S(ok), V(ok), 'o', 'MarkerSize', 4); hold on; grid on;
    plot(ss, polyval(pp, ss), 'r-');
    xlabel('Signal [DN]'); ylabel('Shot-noise variance [DN^2]');
    title('Same, linear axes: variance collapses at saturation');

    figure('Name','Illumination & PRNU','Position',[100 100 1300 450]);
    subplot(1,2,1);
    imagesc(keepMap); axis image off; colorbar; colormap(gca, parula);
    title(sprintf('Light image (open - closed) at %.2f s [DN]', keepT/1e6));
    subplot(1,2,2);
    plot(S(ok), 100*prnu(ok), 'o-', 'MarkerSize', 4); grid on;
    xlabel('Signal [DN]'); ylabel('Spatial std / mean in ROI [%]');
    title('Pixel response non-uniformity (incl. illumination gradient)');

    figure('Name','Dark reference','Position',[130 130 900 420]);
    plot(t, Sc, 'o-', 'MarkerSize', 4); grid on;
    xlabel('Exposure [s]'); ylabel('Closed-shutter mean in ROI [DN]');
    title('Shutter-closed level vs exposure (black-level clamp keeps it ~100 DN)');
end

function [m, vTemp, fullMean] = cubeStats(zipPath, entry, tmp, roiHalf)
    f = imx662_zip_extract(zipPath, entry, tmp);
    c = onCleanup(@() delete(f));
    cube = readCube(f);                               % [ny nx nf] single
    [ny, nx, ~] = size(cube);
    r = round(ny/2) + (-roiHalf+1:roiHalf); q = round(nx/2) + (-roiHalf+1:roiHalf);
    roi = cube(r, q, :);
    m.roi  = mean(roi, 3);
    m.full = mean(cube, 3);
    fullMean = mean(m.full(:));
    d = diff(roi, 1, 3);                              % consecutive-frame differences
    d = d - mean(mean(d, 1), 2);                      % remove frame-to-frame offset/light flicker
    vTemp = mean(var(reshape(d, [], size(d,3)), 0, 1)) / 2;
end

function cube = readCube(fn)
    fid = fopen(fn, 'r', 'ieee-be'); c = onCleanup(@() fclose(fid));
    H = struct(); nb = 0; done = false;
    while ~done
        blk = fread(fid, [80 36], '*char')'; nb = nb + 1;
        for i = 1:36
            card = blk(i,:); key = strtrim(card(1:8));
            if strcmp(key, 'END'), done = true; break; end
            if card(9) == '=', H.(key) = str2double(strtok(card(11:end), '/')); end
        end
    end
    fseek(fid, nb*2880, 'bof');
    cube = fread(fid, H.NAXIS1*H.NAXIS2*H.NAXIS3, 'int16=>single');
    cube = permute(reshape(cube, H.NAXIS1, H.NAXIS2, H.NAXIS3), [2 1 3]);
end
