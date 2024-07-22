% © 2023 Simon Toedtli <s.toedtli@jhu.edu>, all rights reserved
function confirmationTensor = read_confirmation_tensor(dataFolder, fileNrString, ngx, ngy, ngz)
    nElementsTensorBlock = (ngz + 1) * (ngx + 1) * (ngy + 1);
    shapeTensorBlock = [ngz + 1, ngx + 1, ngy + 1];
    fileId = fopen(fullfile(dataFolder, ['Cfield.dble.', fileNrString]));
    rawData = fread(fileId, nElementsTensorBlock, 'double');
    confirmationTensor.Cxx = reshape(rawData, shapeTensorBlock);
    rawData = fread(fileId, nElementsTensorBlock, 'double');
    confirmationTensor.Cyy = reshape(rawData, shapeTensorBlock);
    rawData = fread(fileId, nElementsTensorBlock, 'double');
    confirmationTensor.Czz = reshape(rawData, shapeTensorBlock);
    rawData = fread(fileId, nElementsTensorBlock, 'double');
    confirmationTensor.Cxy = reshape(rawData, shapeTensorBlock);
    rawData = fread(fileId, nElementsTensorBlock, 'double');
    confirmationTensor.Cxz = reshape(rawData, shapeTensorBlock);
    rawData = fread(fileId, nElementsTensorBlock, 'double');
    confirmationTensor.Cyz = reshape(rawData, shapeTensorBlock);
    fclose(fileId);
end
