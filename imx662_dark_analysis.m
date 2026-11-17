%% IMX662 dark-current analysis  (vacuum chamber, 25 C, 2025-07-09)
%
% What is in the dataset
%   11 FITS cubes, one per exposure time (15 us, 0.2, 0.5, 1.0 ... 4.5 s).
%   Each cube = 25 frames of 1920 x 1080, int16 big-endian, 12-bit ADC
%   (saturates at 4095). Analog gain 9, digital gain 1, LED off.
%   Matching *_config_cam.json files hold the acquisition settings.
%
% Key thing to know before reading the plots
%   The sensor/driver clamps the black level to ~100 DN using its optical-
%   black pixels. That clamp subtracts the *average* dark signal, so the
%   median pixel value stays ~100 DN at every exposure. You can't read
%   dark current from the mean level. It shows up in three other places:
%     1) hot pixels: a tail of pixels whose level grows with exposure
%        (dark-signal non-uniformity, DSNU)
%     2) temporal noise: dark electrons bring shot noise, so the
%        frame-to-frame variance rises linearly with exposure time
%     3) the spatial spread (std) of the frame, mostly driven by (1)
%   This script measures all three.
%
% Usage
%   Press Run. A file dialog opens: pick the .zip, or select all .fits files.
%   Results (.mat) are saved in MATLAB's current folder.
%   Needs only base MATLAB (R2018b+). No toolboxes; FITS is read directly.
%   Memory: one cube at a time (~210 MB as single).

% Optional: set dataFile before running to skip the file dialog, e.g.
%   dataFile = '/path/to/IMX662_DarkCurrent_20250709.zip'; imx662_dark_analysis
if ~exist('dataFile','var'), dataFile = ''; end
clearvars -except dataFile; close all; clc;

%% ---------------- User settings ----------------
satLevel  = 4095;     % 12-bit ADC full scale
satGuard  = 4000;     % mean-frame values above this are treated as saturated
hotSigma  = 5;        % hot pixel = slope > median + hotSigma*robust sigma
hotAbsDNs = 10;       % ...and also above this absolute slope [DN/s]
K_eDN     = NaN;      % conversion gain [e-/DN] at this gain setting.
                      % Leave NaN if unknown. Set it (from a photon-
                      % transfer curve) to get results in electrons.

%% ---------------- Pick the data (file dialog) ----------------
% A window opens: select EITHER the .zip you downloaded, OR all the
% .fits files (Ctrl+A / Shift-click to select several at once).
if ~isempty(dataFile)
    [fp, nm, ex] = fileparts(dataFile); fn = {[nm ex]};
else
    [fn, fp] = uigetfile({'*.zip;*.fits;*.fit', 'Dark data (*.zip, *.fits)'; ...
                          '*.*', 'All files'}, ...
                         'Select the IMX662 .zip OR all the .fits files', ...
                         'MultiSelect', 'on');
    if isequal(fn, 0), error('No file selected.'); end
    fn = cellstr(fn);
end

if numel(fn) == 1 && endsWith(lower(fn{1}), '.zip')
    % Unzip into a temporary folder and find every .fits inside it
    unzipDir = fullfile(tempdir, 'imx662_dark_unzipped');
    fprintf('Unzipping %s ...\n', fn{1});
    unzip(fullfile(fp, fn{1}), unzipDir);
    files = dir(fullfile(unzipDir, '**', '*.fits'));
else
    files = dir(fullfile(fp, fn{1}));
    for i = 2:numel(fn), files(i) = dir(fullfile(fp, fn{i})); end
end
files = files(~[files.isdir]);
assert(~isempty(files), 'No .fits files found in the selection.');
nExp = numel(files);

% Read headers first to get exposure times, then sort by exposure
hdrs = cell(nExp,1);
texp = zeros(nExp,1);
for k = 1:nExp
    hdrs{k} = readFitsHeader(fullfile(files(k).folder, files(k).name));
    texp(k) = hdrs{k}.EXPTIME * 1e-6;          % EXPTIME is in microseconds
end
[texp, order] = sort(texp);
files = files(order);  hdrs = hdrs(order);

H  = hdrs{1};
ny = H.NAXIS2; nx = H.NAXIS1; nf = H.NAXIS3;
fprintf('Found %d cubes, %d x %d px, %d frames each\n', nExp, nx, ny, nf);
fprintf('Analog gain %g, digital gain %g, pixel clock %g MHz\n\n', ...
        H.ANAGAIN, H.DIGGAIN, H.PIXCLOCK);

% Per-exposure products
meanMap   = zeros(ny, nx, nExp, 'single');  % average of the 25 frames
varMap    = zeros(ny, nx, nExp, 'single');  % temporal variance per pixel
frameMean = zeros(nf, nExp);                % mean of each frame (clamp stability)
hdrBias   = zeros(nExp,1);  hdrNoise = zeros(nExp,1);

for k = 1:nExp
    fn = fullfile(files(k).folder, files(k).name);
    fprintf('Loading %5.3g s  (%s)\n', texp(k), files(k).name);
    cube = readFitsCube(fn, hdrs{k});              % [ny nx nf] single
    meanMap(:,:,k) = mean(cube, 3);
    varMap(:,:,k)  = var(cube, 0, 3);
    frameMean(:,k) = squeeze(mean(mean(cube,1),2));
    hdrBias(k)  = hdrs{k}.BIAS;
    hdrNoise(k) = hdrs{k}.NOISE;
end
clear cube

%% ---------------- Per-exposure summary statistics ----------------
medLevel   = zeros(nExp,1);  % median of mean frame  -> black level
rawStd     = zeros(nExp,1);  % plain spatial std (hot-pixel dominated)
robStd     = zeros(nExp,1);  % robust spatial std (1.4826*MAD)
medTempVar = zeros(nExp,1);  % median temporal variance
nSat       = zeros(nExp,1);  % pixels saturated in the mean frame
for k = 1:nExp
    m = meanMap(:,:,k);  m = double(m(:));
    medLevel(k)   = median(m);
    rawStd(k)     = std(m);
    robStd(k)     = 1.4826 * median(abs(m - medLevel(k)));
    v = varMap(:,:,k);
    medTempVar(k) = median(double(v(:)));
    nSat(k)       = nnz(m > satGuard);
end

T = table(texp, medLevel, hdrBias, rawStd, robStd, sqrt(medTempVar), ...
          hdrNoise, nSat, 'VariableNames', ...
          {'t_s','median_DN','hdrBIAS','spatialStd_DN','robustStd_DN', ...
           'tempNoise_DN','hdrNOISE','nSaturated'});
disp(T);

%% ---------------- Per-pixel dark-signal slope ----------------
% Remove the clamp: subtract each exposure's median so a typical pixel
% sits at 0. What remains is each pixel's dark signal *relative to the
% array median*. Then fit DN vs t per pixel, skipping saturated points.
y = reshape(meanMap, ny*nx, nExp);                    % [Npix x nExp]
y = double(y) - medLevel.';                            % clamp-corrected
w = double(reshape(meanMap, ny*nx, nExp) < satGuard);  % 1 = usable point
t = texp.';

Sw  = sum(w, 2);
St  = w * texp;          Stt = w * texp.^2;
Sy  = sum(w .* y, 2);    Sty = sum(w .* y .* t, 2);
den = Sw .* Stt - St.^2;
slope     = (Sw .* Sty - St .* Sy) ./ den;   % DN/s
intercept = (Sy - slope .* St) ./ Sw;        % DN at t = 0
slope(Sw < 3) = NaN;   % <3 unsaturated points: not fit (pixel saturates fast)
fastSat = Sw < 3;

slopeMap = reshape(slope, ny, nx);
clear y w

% Robust stats and hot-pixel selection
sv      = slope(~isnan(slope));
sMed    = median(sv);
sSig    = 1.4826 * median(abs(sv - sMed));
hotThr  = max(sMed + hotSigma*sSig, hotAbsDNs);
hotMask = reshape(slope > hotThr | fastSat, ny, nx);

fprintf('\nPer-pixel slope (relative to array median):\n');
fprintf('  median %.3f DN/s, robust sigma %.3f DN/s\n', sMed, sSig);
fprintf('  percentiles 99 / 99.9 / 99.99 %%: %.2f / %.2f / %.2f DN/s\n', ...
        pct(sv, [99 99.9 99.99]));
fprintf('  hot pixels (> %.2f DN/s or saturating fast): %d  (%.4f %% of array)\n', ...
        hotThr, nnz(hotMask), 100*nnz(hotMask)/numel(hotMask));

%% ---------------- Dark shot noise from temporal variance ----------------
% var(t) = readVar + (1/K^2) * D_e * t    (in DN^2, K in e-/DN)
% The slope of median temporal variance vs t gives D_e/K^2 [DN^2/s].
% This is the one handle on the *typical* pixel's dark current, because
% the black-level clamp hides its mean.
pv = polyfit(texp, medTempVar, 1);
fprintf('\nTemporal variance fit: var = %.3f DN^2 + %.3f DN^2/s * t\n', pv(2), pv(1));
fprintf('  read noise (t->0) = %.2f DN rms\n', sqrt(pv(2)));
if ~isnan(K_eDN)
    De = pv(1) * K_eDN^2;
    fprintf('  with K = %.3f e-/DN: typical dark current = %.3f e-/s/px, read noise = %.2f e-\n', ...
            K_eDN, De, sqrt(pv(2))*K_eDN);
else
    fprintf('  (set K_eDN to convert to e-/s)\n');
end

%% ================= FIGURES =================
col = parula(nExp);
lbl = arrayfun(@fmtT, texp, 'UniformOutput', false);

% ---- Fig 1: images (shortest, longest, difference) ----
figure('Name','Dark frames','Position',[50 50 1500 420]);
cl = medLevel(end) + [-5 5];
subplot(1,3,1); imagesc(meanMap(:,:,1), cl); axis image off; colorbar;
title(sprintf('Mean of %d frames, %s', nf, lbl{1}));
subplot(1,3,2); imagesc(meanMap(:,:,end), cl); axis image off; colorbar;
title(sprintf('Mean of %d frames, %s', nf, lbl{end}));
subplot(1,3,3); imagesc(meanMap(:,:,end) - meanMap(:,:,1), [-3 20]); axis image off; colorbar;
title(sprintf('%s minus %s (hot pixels = bright dots)', lbl{end}, lbl{1}));
colormap(gray);

% ---- Fig 2: histograms of mean frame for each exposure ----
figure('Name','Histograms','Position',[80 80 900 500]);
edges = 0:1:satLevel+1;
hold on;
for k = 1:nExp
    c = histcounts(meanMap(:,:,k), edges);
    stairs(edges(1:end-1), max(c,0.5), 'Color', col(k,:), 'LineWidth', 1);
end
set(gca,'YScale','log','XScale','log'); xlim([50 satLevel+1]);
xlabel('Pixel value in mean frame [DN]'); ylabel('Number of pixels');
title('Histogram of mean dark frame: core fixed at ~100 DN, hot-pixel tail grows with t');
legend(lbl, 'Location','northeast'); grid on; box on;

% ---- Fig 3: black level and noise vs exposure ----
figure('Name','Level & noise vs t','Position',[110 110 1300 420]);
subplot(1,3,1);
plot(texp, medLevel, 'o-', texp, hdrBias, 's--', 'LineWidth',1.2); grid on;
xlabel('Exposure [s]'); ylabel('DN');
legend('median of mean frame','header BIAS','Location','best');
title('Black level: flat, so the clamp removes mean dark');
subplot(1,3,2);
plot(texp, rawStd, 'o-', texp, robStd, 's-', 'LineWidth',1.2); grid on;
xlabel('Exposure [s]'); ylabel('Spatial std of mean frame [DN]');
legend('plain std (hot-pixel dominated)','robust std (1.4826*MAD)','Location','northwest');
title('Spatial spread');
subplot(1,3,3);
tt = linspace(0, max(texp), 100);
plot(texp, medTempVar, 'o', tt, polyval(pv, tt), '-', 'LineWidth',1.2); grid on;
xlabel('Exposure [s]'); ylabel('Median temporal variance [DN^2]');
title(sprintf('Dark shot noise: %.2f + %.3f\\cdott  DN^2', pv(2), pv(1)));

% ---- Fig 4: per-pixel slope map and histogram ----
figure('Name','Dark-signal slope','Position',[140 140 1400 480]);
subplot(1,2,1);
imagesc(slopeMap, [-1 5]); axis image off; colorbar; colormap(gca, hot);
title('Per-pixel dark slope [DN/s] relative to array median');
subplot(1,2,2);
se = linspace(min(-3, pct(sv,0.01)), max(pct(sv,99.999), 10), 400);
histogram(sv, se, 'DisplayStyle','stairs','LineWidth',1.2); hold on;
set(gca,'YScale','log'); xline(hotThr,'r--',sprintf('hot threshold %.1f DN/s',hotThr));
xlabel('Slope [DN/s]'); ylabel('Pixels'); grid on;
title(sprintf('Slope histogram: median %.2f, robust \\sigma %.2f DN/s', sMed, sSig));

% ---- Fig 5: response curves of the hottest pixels ----
figure('Name','Hot pixel response','Position',[170 170 1300 450]);
sTmp = slope; sTmp(isnan(sTmp)) = -Inf;
[~, idxHot] = maxk(sTmp, 15);                 % brightest fitted pixels
med1 = find(abs(slope - sMed) < 0.05*sSig, 15); % some typical pixels
Y = reshape(meanMap, ny*nx, nExp);
subplot(1,2,1);
plot(texp, Y(idxHot,:)', 'o-'); grid on; hold on;
yline(satLevel,'k--','saturation');
xlabel('Exposure [s]'); ylabel('Mean-frame value [DN]');
title('15 hottest pixels: linear growth, then saturation at 4095');
subplot(1,2,2);
plot(texp, Y(med1,:)', 'o-'); grid on;
xlabel('Exposure [s]'); ylabel('Mean-frame value [DN]');
title('15 typical pixels: flat at the clamped black level');
clear Y

% ---- Fig 6: clamp stability within each cube ----
figure('Name','Frame-to-frame stability','Position',[200 200 900 450]);
hold on;
for k = 1:nExp, plot(1:nf, frameMean(:,k), '.-', 'Color', col(k,:)); end
xlabel('Frame number'); ylabel('Frame mean [DN]'); grid on; box on;
legend(lbl,'Location','eastoutside');
title('Frame mean vs frame index: black-level clamp drift / settling');

% ---- Fig 7: row / column structure and hot-pixel locations ----
figure('Name','Spatial structure','Position',[230 230 1400 450]);
sm = slopeMap; sm(hotMask) = NaN;
subplot(1,3,1); plot(mean(sm,1,'omitnan')); grid on;
xlabel('Column'); ylabel('Mean slope [DN/s]'); title('Column profile (hot px excluded)');
subplot(1,3,2); plot(mean(sm,2,'omitnan')); grid on;
xlabel('Row'); ylabel('Mean slope [DN/s]'); title('Row profile (hot px excluded)');
subplot(1,3,3);
[r, c] = find(hotMask);
scatter(c, r, 6, slopeMap(hotMask), 'filled'); axis ij image;
xlim([1 nx]); ylim([1 ny]); colorbar; hs = slope(hotMask); hs = hs(~isnan(hs));
if ~isempty(hs), caxis([hotThr max(hotThr+1, pct(hs,95))]); end
title(sprintf('%d hot pixels (color = slope DN/s)', nnz(hotMask)));
xlabel('Column'); ylabel('Row');

% Save the main products
save('imx662_dark_results.mat', 'texp','T','slopeMap','hotMask','hotThr', ...
     'pv','medTempVar','frameMean','sMed','sSig','-v7.3');
fprintf('\nSaved imx662_dark_results.mat\n');

%% ================= Local functions =================
function q = pct(x, p)
% Percentile without the Statistics Toolbox (linear interpolation).
    x = sort(double(x(:))); x = x(~isnan(x)); n = numel(x);
    q = interp1(linspace(0,100,n), x, p);
end

function s = fmtT(t)
    if t < 1e-3, s = sprintf('%.0f \\mus', t*1e6); else, s = sprintf('%.1f s', t); end
end

function H = readFitsHeader(fn)
% Parse the primary FITS header (80-char cards in 2880-byte blocks).
    fid = fopen(fn, 'r'); c = onCleanup(@() fclose(fid));
    H = struct('BSCALE',1,'BZERO',0); nBlocks = 0; done = false;
    while ~done
        blk = fread(fid, [80 36], '*char')'; nBlocks = nBlocks + 1;
        for i = 1:36
            card = blk(i,:);
            key = strtrim(card(1:8));
            if strcmp(key,'END'), done = true; break; end
            if numel(card) < 10 || card(9) ~= '=', continue; end
            val = strtrim(strtok(card(11:end), '/'));
            if any(strcmp(val, {'T','F'}))
                v = strcmp(val,'T');
            elseif ~isempty(val) && val(1) == ''''
                v = strtrim(strrep(val, '''', ''));
            else
                v = str2double(val);
            end
            H.(matlab.lang.makeValidName(key)) = v;
        end
    end
    H.dataOffset = nBlocks * 2880;
end

function cube = readFitsCube(fn, H)
% Read a 3-D FITS image as single, returned as [NAXIS2 NAXIS1 NAXIS3]
% (rows x cols x frames). FITS data is big-endian.
    switch H.BITPIX
        case 8,   p = 'uint8';  case 16, p = 'int16';
        case 32,  p = 'int32';  case -32, p = 'single';
        otherwise, error('Unsupported BITPIX %d', H.BITPIX);
    end
    n = H.NAXIS1 * H.NAXIS2 * H.NAXIS3;
    fid = fopen(fn, 'r', 'ieee-be'); c = onCleanup(@() fclose(fid));
    fseek(fid, H.dataOffset, 'bof');
    cube = fread(fid, n, [p '=>single']);
    cube = reshape(cube, H.NAXIS1, H.NAXIS2, H.NAXIS3);
    cube = permute(cube, [2 1 3]);
    cube = cube * H.BSCALE + H.BZERO;
end
