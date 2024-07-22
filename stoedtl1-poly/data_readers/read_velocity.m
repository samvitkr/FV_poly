% © 2022 Simon Toedtli <s.toedtli@jhu.edu>, all rights reserved
function [uc, vc, wc] = read_velocity(dataFolder, fileNrString, ngx, ngy, ngz)
    % streamwise velocity
    fileId = fopen(fullfile(dataFolder, ['01uc00.dble.', fileNrString]));
    rawData = fread(fileId, 'double');
    uc = reshape(rawData, [ngz+1, ngx+2, ngy+1]);
    fclose(fileId);
    % wall-normal velocity
    fileId = fopen(fullfile(dataFolder, ['02vc00.dble.', fileNrString]));
    rawData = fread(fileId, 'double');
    vc = reshape(rawData, [ngz+1, ngx+1, ngy+2]);
    fclose(fileId);
    % spanwise velocity
    fileId = fopen(fullfile(dataFolder, ['03wc00.dble.', fileNrString]));
    rawData = fread(fileId, 'double');
    wc = reshape(rawData, [ngz+2, ngx+1, ngy+1]);
    fclose(fileId);
end

