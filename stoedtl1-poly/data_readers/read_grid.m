% © 2022 Simon Toedtli <s.toedtli@jhu.edu>, all rights reserved
function [xGridPoint, yGridPoint, xCenter, yCenter] = read_grid(dataFolder, ngx, ngy)
    fileId = fopen(fullfile(dataFolder, 'grid.data'));
    inputBuffer = fread(fileId, 'double');
    xMeshGrid = reshape(inputBuffer(1 : ngx*ngy), [ngx, ngy]);
    yMeshGrid= reshape(inputBuffer(ngx*ngy+1 : 2*ngx*ngy), [ngx, ngy]);
    xGridPoint = squeeze(xMeshGrid(:, 1));  % same x coordinate vector at all y, pick first one
    yGridPoint = squeeze(yMeshGrid(1, :));  % same y coordinate vector at all x, pick first one
    xCenter = 0.5 * (xGridPoint(1 : ngx-1) + xGridPoint(2 : ngx));
    yCenter = 0.5 * (yGridPoint(1 : ngy-1) + yGridPoint(2 : ngy));
    fclose(fileId);
end