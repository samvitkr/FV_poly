clear;  clc;  close all;
addpath('../data_readers/');

runName = 'test_wi_6-67_1';
fileNr = 10;
ngx = 193;
ngy = 225;
ngz = 161;
Lz = 2.0 * pi;
dataFolder = fullfile(getenv('HOME'), 'research_data', 'polymer');

% read grid and generate Fourier transform object
runFolder = fullfile(dataFolder, runName);
[xGridPointsDns, ~, ~, ~] = read_grid(runFolder, ngx, ngy);
deltaX = xGridPointsDns(2) - xGridPointsDns(1);
deltaZ = Lz / (ngz - 1);
ft = SpatialFourierTransform();

% read confirmation tensor
fileNrString = num2str(fileNr, '%07d');
confirmationTensor = read_confirmation_tensor(runFolder, fileNrString, ngx, ngy, ngz);
confirmationTensor = remove_confirmation_tensor_ghost_points(confirmationTensor);

% transform to Fourier domain and back
[confirmationTensorFourier, kx, kz] = ft.transform_confirmation_tensor_to_fourier(confirmationTensor, deltaX, deltaZ);
confirmationTensorAlt = ft.transform_confirmation_tensor_to_physical(confirmationTensorFourier);

% make sure the twice-transformed fields are identical up to machine precision
diffXX = confirmationTensor.Cxx - confirmationTensorAlt.Cxx;
diffYY = confirmationTensor.Cyy - confirmationTensorAlt.Cyy;
diffZZ = confirmationTensor.Czz - confirmationTensorAlt.Czz;
diffXY = confirmationTensor.Cxy - confirmationTensorAlt.Cxy;
diffXZ = confirmationTensor.Cxz - confirmationTensorAlt.Cxz;
diffYZ = confirmationTensor.Cyz - confirmationTensorAlt.Cyz;
disp(['Maximum value Cxx = ', num2str(get_global_max(confirmationTensor.Cxx)), ', maximum difference = ', ...
    num2str(get_global_max(diffXX))]);
disp(['Maximum value Cyy = ', num2str(get_global_max(confirmationTensor.Cyy)), ', maximum difference = ', ...
    num2str(get_global_max(diffYY))]);
disp(['Maximum value Czz = ', num2str(get_global_max(confirmationTensor.Czz)), ', maximum difference = ', ...
    num2str(get_global_max(diffZZ))]);
disp(['Maximum value Cxy = ', num2str(get_global_max(confirmationTensor.Cxy)), ', maximum difference = ', ...
    num2str(get_global_max(diffXY))]);
disp(['Maximum value Cxz = ', num2str(get_global_max(confirmationTensor.Cxz)), ', maximum difference = ', ...
    num2str(get_global_max(diffXZ))]);
disp(['Maximum value Cyz = ', num2str(get_global_max(confirmationTensor.Cyz)), ', maximum difference = ', ...
    num2str(get_global_max(diffYZ))]);


function qMax = get_global_max(q)
    qMax = max(abs(q), [], 'all');
end
