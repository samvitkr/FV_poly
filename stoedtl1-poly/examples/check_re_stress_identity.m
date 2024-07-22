% © 2023 Simon Toedtli <s.toedtli@jhu.edu>, all rights reserved
clear;  clc;  close all;

% automatically add code modules to matlab path
functionCallStack = dbstack;
[scriptFolder, ~, ~] = fileparts(which(functionCallStack(1).file));
[sourceFolder, ~, ~] = fileparts(scriptFolder);
addpath(fullfile(sourceFolder, 'data_readers'), fullfile(sourceFolder, 'utilities'));

% user input: specify run parameters. This assumes your data is stored in a directory with full path dataFolder/runName
runName = 'canonical_4pi_2pi_5';
fileNr = 400;
ngx = 193;  % ng{i}: grid size of restart files in coordinate i, including periodic points in x, z.
ngy = 225;  %        ghost points do not have to be included here, they are accounted for in the data reading routines
ngz = 161;
nCheb = 151;  % number of Chebyshev collocation points in y (has to be determined a priori)
Lz = 2.0 * pi;
dataFolder = fullfile(getenv('HOME'), 'research_data');

% read grid and generate operators
runFolder = fullfile(dataFolder, runName);
[xGridPointsDns, yGridPointsDns, ~, ~] = read_grid(runFolder, ngx, ngy);
deltaX = xGridPointsDns(2) - xGridPointsDns(1);
deltaZ = Lz / (ngz - 1);
zGridPointsDns = 0 : deltaZ : Lz;
[yCheb, ~] = chebdif(nCheb, 1);
interpY = InterpolationOperatorY(yGridPointsDns, yCheb);
ft = SpatialFourierTransform();

% read DNS data: cartesian components in phyiscal domain, represented on staggered grid
fileNrString = num2str(fileNr, '%07d');
[u, v, w] = read_velocity(runFolder, fileNrString, ngx, ngy, ngz);  % these arrays still contain ghost points
[u, v, w] = remove_velocity_ghost_points(u, v, w);

% collocate data spectrally in x, z
[uFourier, kx, kz] = ft.transform_to_fourier(u, deltaX, deltaZ);
vFourier = ft.transform_to_fourier(v);  % no need to pass sampling rate, wavenumber object is needed only once
wFourier = ft.transform_to_fourier(w);
uFourier = interpolate_fourier_in_z(uFourier, kz, -0.5*deltaZ);
vFourier = interpolate_fourier_in_x(vFourier, kx, -0.5*deltaX);
vFourier = interpolate_fourier_in_z(vFourier, kz, -0.5*deltaZ);
wFourier = interpolate_fourier_in_x(wFourier, kx, -0.5*deltaX);

% collocate on Chebyshev grid in y using spline interpolation
uFourier = interpY.interpolate_u_to_chebyshev_grid(uFourier);
vFourier = interpY.interpolate_v_to_chebyshev_grid(vFourier);
wFourier = interpY.interpolate_w_to_chebyshev_grid(wFourier);

% the velocity fields were divergence-free w.r.t. the finite volume numerics, but not necessarily w.r.t the spectral
% numerics. We therefore project onto a divergence free subspace (note: the projected fields satisfy no-through at the
% wall, but they may have non-zero slip).
[uFourier, vFourier, wFourier] = make_field_divergence_free(uFourier, vFourier, wFourier, kx, kz, nCheb);

% compute vorticity components
dx = FirstDerivativeXFourier(kx.get_radial_frequency_vector());
dz = FirstDerivativeZFourier(kz.get_radial_frequency_vector());
dy = FirstDerivativeYChebyshev(nCheb);
omegaYFourier = dz.compute_derivative(uFourier) - dx.compute_derivative(wFourier);
omegaZFourier = dx.compute_derivative(vFourier) - dy.compute_derivative(uFourier);

% map back to physical space, compute averages and flow identity residual
u = ft.transform_to_physical(uFourier);
v = ft.transform_to_physical(vFourier);
w = ft.transform_to_physical(wFourier);
omegaY = ft.transform_to_physical(omegaYFourier);
omegaZ = ft.transform_to_physical(omegaZFourier);
meanReynoldsStress = dy.compute_derivative(squeeze(mean(u .* v, [1 2])));
meanRhs1 = squeeze(mean(w .* omegaY, [1 2]));
meanRhs2 = squeeze(mean(v .* omegaZ, [1 2]));
residual = abs(meanReynoldsStress - (meanRhs1 - meanRhs2));

% plot terms and residual
figure();  % individual terms
plot(yCheb, meanReynoldsStress, 'LineWidth', 2.0, 'DisplayName', '$\mathrm{d}_y \overline{uv}$');
hold on;
plot(yCheb, meanRhs1, 'LineWidth', 2.0, 'DisplayName', '$\overline{w \omega_y}$');
plot(yCheb, -meanRhs2, 'LineWidth', 2.0, 'DisplayName', '$-\overline{v \omega_z}$');
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('$y/h$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$q \, (h / U_b^2)$', 'Interpreter', 'latex', 'FontSize', 14);
title('Quantities entering flow identity')

rgbOrange = [0.850, 0.325, 0.098];
figure();  % residual
plot(yCheb, residual, 'Color', rgbOrange, 'LineWidth', 2.0);
xlabel('$y/h$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$|\mathrm{d}_y \overline{uv} - \overline{w \omega_y} + \overline{v \omega_z}| \, (h / U_b^2)$', ...
    'Interpreter', 'latex', 'FontSize', 14);
title('Residual flow identity')

