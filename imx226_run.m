% imx226_run.m -- full IMX226 characterization (Carissma 2026), same lab method as the IMX662.
%   Inputs (in Downloads):
%     Box_data/IMX226_reduced/IMX226_Po210, .../IMX226_Fe55   reduced event data (imx226_box_reducer.js)
%   Stages: 1 Po-210, 6 Fe-55, 2 readout noise, 3 optical linearity, 4 dark current, 5 tables, 7 beamline
%     Box_data/IMX226_reduced/IMX226_readnoise                 per-pixel readout noise (bias maps)
%     Box_data/IMX226_reduced/IMX226_linearity                 per-frame optical linearity statistics
%     IMX226_DarkCurrent_20250604.zip                          dark-current cubes (full zip from Box)
%   Output: Downloads/IMX226_results (PNG + FIG + CSV + .mat)
close all; clear functions;
c = imx662_cfg();
RD  = fullfile(c.D, 'Box_data', 'IMX226_reduced');
OUT = fullfile(c.D, 'IMX226_results'); if ~exist(OUT, 'dir'), mkdir(OUT); end
diary(fullfile(OUT, 'imx226_log.txt')); diary on
if ~exist('stages', 'var'), stages = [1 6 2 3 4 5]; end
fprintf('\n######## IMX226 run %s, stages %s ########\n', datestr(now), mat2str(stages));

if any(stages == 1)   % X-ray response
    Rpo = imx226_xray_from_reduced(fullfile(RD, 'IMX226_Po210'), 'Po210');
    save(fullfile(OUT, 'imx226_po210.mat'), 'Rpo', '-v7.3');
    imx662_plot_xray(Rpo); imx662_save_figs('IMX226_Po210', OUT);
end
if any(stages == 6)   % Fe-55 X-ray response
    if ~exist('Rpo', 'var'), load(fullfile(OUT, 'imx226_po210.mat'), 'Rpo'); end
    Rfe = imx226_xray_from_reduced(fullfile(RD, 'IMX226_Fe55'), 'Fe55', Rpo.gain);
    save(fullfile(OUT, 'imx226_fe55.mat'), 'Rfe', '-v7.3');
    imx662_plot_xray(Rfe); imx662_save_figs('IMX226_Fe55', OUT);
end
if any(stages == 2)   % readout noise
    if ~exist('Rpo', 'var'), load(fullfile(OUT, 'imx226_po210.mat'), 'Rpo'); end
    RN = imx226_readnoise(fullfile(RD, 'IMX226_readnoise'), Rpo.gain);
    save(fullfile(OUT, 'imx226_readnoise.mat'), 'RN'); imx662_save_figs('IMX226_ReadoutNoise', OUT);
end
if any(stages == 3)   % optical linearity
    if ~exist('Rpo', 'var'), load(fullfile(OUT, 'imx226_po210.mat'), 'Rpo'); end
    OL = imx226_optical_linearity(fullfile(RD, 'IMX226_linearity'), Rpo.gain);
    save(fullfile(OUT, 'imx226_optical_linearity.mat'), 'OL'); imx662_save_figs('IMX226_OpticalLinearity', OUT);
end
if any(stages == 4)   % dark current
    if ~exist('Rpo', 'var'), load(fullfile(OUT, 'imx226_po210.mat'), 'Rpo'); end
    DC = imx226_dark_current(Rpo.gain);
    save(fullfile(OUT, 'imx226_dark.mat'), 'DC'); imx662_save_figs('IMX226_Dark', OUT);
end
if any(stages == 5)   % X-ray linearity, resolution budget and tables
    imx226_tables(OUT);
    imx662_save_figs('IMX226_Tables', OUT);
end
if any(stages == 7)   % polarimetry beamline background cubes
    if ~exist('Rpo', 'var'), load(fullfile(OUT, 'imx226_po210.mat'), 'Rpo'); end
    BL = imx226_beamline(fullfile(RD, 'IMX226_beamline'), Rpo.gain);
    save(fullfile(OUT, 'imx226_beamline.mat'), 'BL', '-v7.3'); imx662_save_figs('IMX226_Beamline', OUT);
end
fprintf('######## IMX226 finished %s ########\n', datestr(now));
diary off
disp('IMX226 RUN DONE')
