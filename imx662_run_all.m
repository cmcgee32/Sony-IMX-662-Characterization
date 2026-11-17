%% imx662_run_all.m -- full IMX662 characterization in the lab's standard format.
% Produces every IMX662 plot/table type of Lupo (2025) Ch. 4 for the four datasets,
% plus extra plots, each with a plain-language caption.
%   stage 1 = dark current        (IMX662_DarkCurrent_20250709.zip)
%   stage 2 = optical linearity   (IMX662_Linearity_20250624.zip)
%   stage 3 = Fe-55 X-ray response, hot pixels, light curves, heat maps, gain variance
%   stage 4 = Po-210 X-ray response (C, O, F lines)
%   stage 5 = readout noise (noise maps from both X-ray runs)
%   stage 6 = X-ray linearity, energy-resolution budget, all tables
%   stage 7 = only re-draw the Fe-55 / Po-210 figures from saved .mat results
% Usage:  stages = [1 2 3 4 5 6]; imx662_run_all     (default: all)
% Outputs: Downloads/IMX662_results/*.png, *.fig, *.csv, *.mat and imx662_log.txt

if ~exist('stages', 'var'), stages = 1:6; end
cfg = imx662_cfg(); OUT = cfg.OUT;
diary(fullfile(OUT, 'imx662_log.txt')); diary on;
fprintf('\n######## IMX662 run started %s, stages %s ########\n', datestr(now), mat2str(stages));
gainFile = fullfile(OUT, 'fe55_results.mat');
getGain = @() cfg.gainRef;
if exist(gainFile, 'file'), tmpG = load(gainFile, 'gain'); getGain = @() tmpG.gain; end

if ismember(3, stages)
    fprintf('\n################ 3. Fe-55 ################\n');
    Fe = imx662_xray_events('Fe55', cfg.gainRef);
    save(gainFile, '-struct', 'Fe', '-v7.3'); tmpG.gain = Fe.gain; getGain = @() tmpG.gain;
    imx662_plot_xray(Fe); imx662_save_figs('3_Fe55', OUT);
end
if ismember(4, stages)
    fprintf('\n################ 4. Po-210 ################\n');
    Po = imx662_xray_events('Po210', getGain());
    save(fullfile(OUT, 'po210_results.mat'), '-struct', 'Po', '-v7.3');
    imx662_plot_xray(Po); imx662_save_figs('4_Po210', OUT);
end
if ismember(7, stages)      % re-draw the X-ray figures from saved results (no re-processing)
    Fe = load(gainFile); imx662_plot_xray(Fe); imx662_save_figs('3_Fe55', OUT);
    Po = load(fullfile(OUT, 'po210_results.mat')); imx662_plot_xray(Po); imx662_save_figs('4_Po210', OUT);
end
if ismember(1, stages)
    fprintf('\n################ 1. DARK CURRENT ################\n');
    DC = imx662_dark_current(getGain());
    save(fullfile(OUT, 'dark_results.mat'), '-struct', 'DC', '-v7.3');
    imx662_save_figs('1_Dark', OUT);
end
if ismember(2, stages)
    fprintf('\n################ 2. OPTICAL LINEARITY ################\n');
    OL = imx662_optical_linearity(getGain());
    save(fullfile(OUT, 'optical_linearity_results.mat'), '-struct', 'OL');
    imx662_save_figs('2_OpticalLinearity', OUT);
end
if ismember(5, stages)
    fprintf('\n################ 5. READOUT NOISE ################\n');
    RN = imx662_readnoise(getGain());
    save(fullfile(OUT, 'readnoise_results.mat'), '-struct', 'RN', '-v7.3');
    imx662_save_figs('5_ReadoutNoise', OUT);
end
if ismember(6, stages)
    fprintf('\n################ 6. X-RAY LINEARITY, RESOLUTION BUDGET, TABLES ################\n');
    Fe = load(gainFile); Po = load(fullfile(OUT, 'po210_results.mat'));
    RN = load(fullfile(OUT, 'readnoise_results.mat')); OL = load(fullfile(OUT, 'optical_linearity_results.mat'));
    DC = load(fullfile(OUT, 'dark_results.mat'));
    XL = imx662_xray_linearity(Fe, Po, RN);
    save(fullfile(OUT, 'xray_linearity_results.mat'), '-struct', 'XL');
    imx662_save_figs('6_XrayLinearity', OUT);
    imx662_tables(Fe, Po, RN, XL, OL, DC);
end
fprintf('\n######## finished %s ########\n', datestr(now));
diary off;
disp('ALL REQUESTED STAGES DONE');
