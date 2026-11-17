function names = imx662_zip_list(zipPath, pattern)
% Names of zip entries whose name matches the regular expression 'pattern'.
    zf = java.util.zip.ZipFile(zipPath);
    c = onCleanup(@() zf.close());
    en = zf.entries(); names = {};
    while en.hasMoreElements()
        e = en.nextElement(); n = char(e.getName());
        if ~e.isDirectory() && ~isempty(regexp(n, pattern, 'once'))
            names{end+1,1} = n; %#ok<AGROW>
        end
    end
end
