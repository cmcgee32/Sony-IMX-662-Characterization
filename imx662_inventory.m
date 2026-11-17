% imx662_inventory.m  -- list what is inside each IMX662 zip without unzipping.
% Writes imx662_inventory.txt next to the zips.
d = '/Users/carissmamcgee/Downloads';
zips = {'IMX662_DarkCurrent_20250709.zip','IMX662_Linearity_20250624.zip', ...
        'IMX662_Fe55_SpectrumData_20250625.zip','IMX662_Po210_SpectrumData_20250626.zip'};
out = fopen(fullfile(d,'imx662_inventory.txt'),'w');
for z = 1:numel(zips)
    zp = fullfile(d, zips{z});
    fprintf(out, '\n==================== %s ====================\n', zips{z});
    zf = java.util.zip.ZipFile(zp);
    en = zf.entries(); n = 0; shownJson = 0; shownFits = 0; shownOther = 0;
    while en.hasMoreElements()
        e = en.nextElement(); nm = char(e.getName()); sz = e.getSize(); n = n + 1;
        fprintf(out, '%12d  %s\n', sz, nm);
        low = lower(nm);
        isTxt = endsWith(low,'.json') || endsWith(low,'.txt') || endsWith(low,'.csv') || endsWith(low,'.md');
        if ~e.isDirectory() && isTxt && sz < 200000 && shownJson < 4
            fprintf(out, '----- contents of %s -----\n%s\n-----\n', nm, readEntryText(zf, e, 4000));
            shownJson = shownJson + 1;
        elseif endsWith(low,'.fits') && shownFits < 3
            fprintf(out, '----- FITS header of %s -----\n%s-----\n', nm, readEntryText(zf, e, 2880*2));
            shownFits = shownFits + 1;
        elseif ~e.isDirectory() && ~endsWith(low,'.fits') && ~isTxt && shownOther < 3 && sz < 5e6
            shownOther = shownOther + 1;   % just listed
        end
    end
    fprintf(out, 'TOTAL ENTRIES: %d\n', n);
    zf.close();
end
fclose(out);
disp('INVENTORY DONE');

function txt = readEntryText(zf, e, nmax)
    is = zf.getInputStream(e);
    buf = zeros(1, nmax, 'uint8'); k = 0;
    while k < nmax
        b = is.read();
        if b < 0, break; end
        k = k + 1; buf(k) = b;
    end
    is.close();
    txt = char(buf(1:k));
    if all(txt(1:min(end,6)) == 'SIMPLE') % FITS header: split into 80-char cards
        nc = floor(numel(txt)/80); cards = reshape(txt(1:nc*80), 80, nc)';
        iend = find(strncmp(cellstr(cards), 'END', 3) & all(cards(:,4:end) == ' ', 2), 1);
        if ~isempty(iend), cards = cards(1:iend,:); end
        keep = ~all(cards == ' ', 2);
        txt = sprintf('%s\n', string(cellstr(cards(keep,:))));
    end
end
