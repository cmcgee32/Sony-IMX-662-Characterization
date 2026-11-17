% imx662_probe.m -- extract one example of each event-mode file and describe it.
d = '/Users/carissmamcgee/Downloads';
tmp = fullfile(tempdir, 'imx662_probe'); if exist(tmp,'dir'), rmdir(tmp,'s'); end; mkdir(tmp);
out = fopen(fullfile(d,'imx662_probe.txt'),'w');
items = {
 'IMX662_Fe55_SpectrumData_20250625.zip','20250625/20250625_221311IMX662_VacuumChamber_Fe55__1_evts.fits'
 'IMX662_Fe55_SpectrumData_20250625.zip','20250625/20250625_221311IMX662_VacuumChamber_Fe55__1_grades.fits'
 'IMX662_Fe55_SpectrumData_20250625.zip','20250625/20250625_221311IMX662_VacuumChamber_Fe55__1_biasmap.fits'
 'IMX662_Fe55_SpectrumData_20250625.zip','20250625/hot_pixels.fits'
 'IMX662_Po210_SpectrumData_20250626.zip','20250626/20250626_221924IMX662_VacuumChamber_Po_406_evts.fits'
 'IMX662_Po210_SpectrumData_20250626.zip','20250626/20250626_221924IMX662_VacuumChamber_Po_405_grades.fits'
 };
for i = 1:size(items,1)
    zp = fullfile(d, items{i,1}); en = items{i,2};
    try
        f = extractOne(zp, en, tmp);
        fprintf(out, '\n######## %s :: %s\n', items{i,1}, en);
        info = fitsinfo(f);
        fprintf(out, 'Contents: %s\n', strjoin(info.Contents, ' | '));
        % primary
        kw = info.PrimaryData.Keywords;
        for r = 1:size(kw,1), fprintf(out, 'P  %-8s = %s\n', kw{r,1}, val2str(kw{r,2})); end
        % image extensions
        if isfield(info,'Image')
            for e = 1:numel(info.Image)
                kw = info.Image(e).Keywords; fprintf(out, '--- Image ext %d size %s\n', e, mat2str(info.Image(e).Size));
                for r = 1:size(kw,1), fprintf(out, 'I%d %-8s = %s\n', e, kw{r,1}, val2str(kw{r,2})); end
                im = fitsread(f, 'image', e);
                v = double(im(:)); v = v(isfinite(v));
                fprintf(out, 'I%d stats: min %g  median %g  mean %g  max %g  nnz %d of %d\n', e, min(v), median(v), mean(v), max(v), nnz(v), numel(v));
            end
        end
        if isfield(info,'BinaryTable')
            for e = 1:numel(info.BinaryTable)
                bt = info.BinaryTable(e); kw = bt.Keywords;
                fprintf(out, '--- BinTable ext %d rows %d\n', e, bt.Rows);
                for r = 1:size(kw,1), fprintf(out, 'B%d %-8s = %s\n', e, kw{r,1}, val2str(kw{r,2})); end
                t = fitsread(f, 'binarytable', e, 'TableRows', 1:min(8, bt.Rows));
                for c = 1:numel(t)
                    x = t{c}; fprintf(out, '  col %d (%s) size %s first rows: %s\n', c, colname(kw,c), mat2str(size(x)), mat2str(double(x(1:min(end,24))),5));
                end
            end
        end
        delete(f);
    catch ME
        fprintf(out, 'ERROR on %s: %s\n', en, ME.message);
    end
end
fclose(out); rmdir(tmp,'s'); disp('PROBE DONE');

function s = val2str(v)
    if ischar(v), s = v; elseif islogical(v), s = mat2str(v); else, s = num2str(v, 10); end
end
function n = colname(kw, c)
    k = find(strcmp(kw(:,1), sprintf('TTYPE%d', c)), 1); if isempty(k), n = '?'; else, n = kw{k,2}; end
end
function outFile = extractOne(zipPath, entryName, outDir)
    [~, nm, ex] = fileparts(entryName); outFile = fullfile(outDir, [nm ex]);
    st = system(sprintf('unzip -p "%s" "%s" > "%s"', zipPath, entryName, outFile));
    assert(st == 0, 'unzip failed for %s', entryName);
end
