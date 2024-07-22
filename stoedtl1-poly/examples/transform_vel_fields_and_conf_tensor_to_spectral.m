% © 2023 Simon Toedtli <s.toedtli@jhu.edu>, all rights reserved
clear;  clc;  close all;
addpath 'C:\Users\samvi\Dropbox\SimonsProject\Finite_Vol_data\dmsuite'

% automatically add code modules to matlab path
functionCallStack = dbstack;
[scriptFolder, ~, ~] = fileparts(which(functionCallStack(1).file));
[sourceFolder, ~, ~] = fileparts(scriptFolder);
addpath(fullfile(sourceFolder, 'data_readers'), fullfile(sourceFolder, 'utilities'));

% user input: specify run parameters. This assumes your data is stored in a directory with full path dataFolder/runName
runName = 'test_wi_6-67_1';
fileNr = 79000;
ngx = 129;  % ng{i}: grid size of restart files in coordinate i, including periodic points in x, z.
ngy = 201;  %        ghost points do not have to be included here, they are accounted for in the data reading routines
ngz = 129;
nCheb = 151;  % number of Chebyshev collocation points in y (has to be determined a priori)
Lz = 2.0 * pi;
dataFolder = fullfile(getenv('HOME'), 'research_data', 'polymer');

% read grid and generate operators
%runFolder = fullfile(dataFolder, runName);
runFolder= 'C:\Users\samvi\Dropbox\SimonsProject\Finite_Vol_data\Visco_data_sample';

[xGridPointsDns, yGridPointsDns, ~, ~] = read_grid(runFolder, ngx, ngy);
deltaX = xGridPointsDns(2) - xGridPointsDns(1);
deltaZ = Lz / (ngz - 1);
[yCheb, ~] = chebdif(nCheb, 1);
interpY = InterpolationOperatorY(yGridPointsDns, yCheb);
ft = SpatialFourierTransform();

% read DNS data: Cartesian velocities and confirmation tensor in physical domain, represented on staggered grid
fileNrString = num2str(fileNr, '%07d');
[u, v, w] = read_velocity(runFolder, fileNrString, ngx, ngy, ngz);  % these arrays still contain ghost points
confTensor = read_confirmation_tensor(runFolder, fileNrString, ngx, ngy, ngz);
[u, v, w] = remove_velocity_ghost_points(u, v, w);
confTensor = remove_confirmation_tensor_ghost_points(confTensor);

% tranform to Fourier domain
[uFourier, kx, kz] = ft.transform_to_fourier(u, deltaX, deltaZ);
vFourier = ft.transform_to_fourier(v);  % no need to pass sampling rate, wavenumber object is needed only once
wFourier = ft.transform_to_fourier(w);
confTensorFourier = ft.transform_confirmation_tensor_to_fourier(confTensor);

% collocate data spectrally in x, z
shiftX = -0.5 * deltaX;
shiftZ = -0.5 * deltaZ;
uFourier = interpolate_fourier_in_z(uFourier, kz, shiftZ);
vFourier = interpolate_fourier_in_x(vFourier, kx, shiftX);
vFourier = interpolate_fourier_in_z(vFourier, kz, shiftZ);
wFourier = interpolate_fourier_in_x(wFourier, kx, shiftX);
confTensorFourier = interpolate_confirmation_tensor_fourier_in_x(confTensorFourier, kx, shiftX);
confTensorFourier = interpolate_confirmation_tensor_fourier_in_z(confTensorFourier, kz, shiftZ);

% collocate on Chebyshev grid in y using spline interpolation
uFourier = interpY.interpolate_u_to_chebyshev_grid(uFourier);
vFourier = interpY.interpolate_v_to_chebyshev_grid(vFourier);
wFourier = interpY.interpolate_w_to_chebyshev_grid(wFourier);
confTensorFourier = interpY.interpolate_confirmation_tensor_to_chebyshev_grid(confTensorFourier);
confTensorPhysical=ft.transform_confirmation_tensor_to_physical(confTensorFourier);
% the velocity fields were divergence-free w.r.t. the finite volume numerics, but not necessarily w.r.t the spectral
% numerics. We therefore project onto a divergence free subspace (note: the projected fields satisfy no-through at the
% wall, but they may have non-zero slip).
[uFourier, vFourier, wFourier] = make_field_divergence_free(uFourier, vFourier, wFourier, kx, kz, nCheb);
