function outFile = imx662_zip_extract(zipPath, entryName, outDir)
% Extract ONE entry of a zip to outDir (macOS/Linux 'unzip', fast; no full unzip).
% The caller deletes outFile when done, so only one file uses disk at a time.
    [~, nm, ex] = fileparts(entryName);
    outFile = fullfile(outDir, [nm ex]);
    st = system(sprintf('unzip -p "%s" "%s" > "%s"', zipPath, entryName, outFile));
    if st ~= 0
        if exist(outFile, 'file'), delete(outFile); end
        error('Could not extract %s (disk full?)', entryName);
    end
end
