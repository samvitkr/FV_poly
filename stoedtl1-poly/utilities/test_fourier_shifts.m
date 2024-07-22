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

% read velocity field and generate a mock-up confirmation tensor to test the new shift routines
fileNrString = num2str(fileNr, '%07d');
[u, v, w] = read_velocity(runFolder, fileNrString, ngx, ngy, ngz);
[u, v, w] = remove_velocity_ghost_points(u, v, w);
confirmationTensor.Cxx = v;
confirmationTensor.Cyy = v;
confirmationTensor.Czz = v;
confirmationTensor.Cxy = v;
confirmationTensor.Cxz = v;
confirmationTensor.Cyz = v;

% transform to Fourier domain
[vFourier, kx, kz] = ft.transform_to_fourier(v, deltaX, deltaZ);
confirmationTensorFourier = ft.transform_confirmation_tensor_to_fourier(confirmationTensor);

% apply shift in x and check error
shiftX = -0.5 * deltaX;
vFourier = interpolate_fourier_in_x(vFourier, kx, shiftX);
confirmationTensorFourier = interpolate_confirmation_tensor_fourier_in_x(confirmationTensorFourier, kx, shiftX);
disp('Difference after shift in x');
disp('---------------------------');
diffXX = confirmationTensorFourier.Cxx - vFourier;
display_global_max(diffXX, 'Cxx');
diffYY = confirmationTensorFourier.Cyy - vFourier;
display_global_max(diffYY, 'Cyy');
diffZZ = confirmationTensorFourier.Czz - vFourier;
display_global_max(diffZZ, 'Czz');
diffXY = confirmationTensorFourier.Cxy - vFourier;
display_global_max(diffXY, 'Cxy');
diffXZ = confirmationTensorFourier.Cxz - vFourier;
display_global_max(diffXZ, 'Cxz');
diffYZ = confirmationTensorFourier.Cyz - vFourier;
display_global_max(diffYZ, 'Cyz');

% apply shift in z and check error
shiftZ = -0.5 * deltaZ;
vFourier = interpolate_fourier_in_z(vFourier, kz, shiftZ);
confirmationTensorFourier = interpolate_confirmation_tensor_fourier_in_z(confirmationTensorFourier, kz, shiftZ);
disp('Difference after shift in z');
disp('---------------------------');
diffXX = confirmationTensorFourier.Cxx - vFourier;
display_global_max(diffXX, 'Cxx');
diffYY = confirmationTensorFourier.Cyy - vFourier;
display_global_max(diffYY, 'Cyy');
diffZZ = confirmationTensorFourier.Czz - vFourier;
display_global_max(diffZZ, 'Czz');
diffXY = confirmationTensorFourier.Cxy - vFourier;
display_global_max(diffXY, 'Cxy');
diffXZ = confirmationTensorFourier.Cxz - vFourier;
display_global_max(diffXZ, 'Cxz');
diffYZ = confirmationTensorFourier.Cyz - vFourier;
display_global_max(diffYZ, 'Cyz');


function qMax = get_global_max(q)
    qMax = max(abs(q), [], 'all');
end

function display_global_max(q, qStr)
    disp(['Maximum deviation in ', qStr, ' = ', num2str(get_global_max(q))]);
end
