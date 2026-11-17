function cube = imx662_read_cube(fn)
% Read a 3-D FITS image cube (int16, big-endian) as single [rows x cols x frames].
    fid = fopen(fn, 'r', 'ieee-be'); c = onCleanup(@() fclose(fid));
    H = struct('BSCALE', 1, 'BZERO', 0); nb = 0; done = false;
    while ~done
        blk = fread(fid, [80 36], '*char')'; nb = nb + 1;
        for i = 1:36
            card = blk(i,:); key = strtrim(card(1:8));
            if strcmp(key, 'END'), done = true; break; end
            if card(9) == '=' && isvarname(key), H.(key) = str2double(strtok(card(11:end), '/')); end
        end
    end
    fseek(fid, nb*2880, 'bof');
    cube = fread(fid, H.NAXIS1*H.NAXIS2*H.NAXIS3, 'int16=>single');
    cube = permute(reshape(cube, H.NAXIS1, H.NAXIS2, H.NAXIS3), [2 1 3]);
    if H.BSCALE ~= 1 || H.BZERO ~= 0, cube = cube * H.BSCALE + H.BZERO; end
end
