% go.m -- run the full IMX662 analysis (type "go" in MATLAB).
close all; clear functions;
OUTd = '/Users/carissmamcgee/Downloads/IMX662_results';
old = [dir(fullfile(OUTd, '1_dark_*')); dir(fullfile(OUTd, '1_Dark_*'))];
for k = 1:numel(old), delete(fullfile(OUTd, old(k).name)); end   % remove earlier drafts
stages = [3 4 1 2 5 6];
imx662_run_all
