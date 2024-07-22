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
reBulk = 2800;
reTau = 180;
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

% compute vorticity and requried gradients
dx = FirstDerivativeXFourier(kx.get_radial_frequency_vector());
dz = FirstDerivativeZFourier(kz.get_radial_frequency_vector());
dy = FirstDerivativeYChebyshev(nCheb);
omegaXFourier = dy.compute_derivative(wFourier) - dz.compute_derivative(vFourier);
omegaYFourier = dz.compute_derivative(uFourier) - dx.compute_derivative(wFourier);
omegaZFourier = dx.compute_derivative(vFourier) - dy.compute_derivative(uFourier);
dzOmegaYFourier = dz.compute_derivative(omegaYFourier);

% map back to physical space
u = ft.transform_to_physical(uFourier);
v = ft.transform_to_physical(vFourier);
w = ft.transform_to_physical(wFourier);
omegaX = ft.transform_to_physical(omegaXFourier);
omegaY = ft.transform_to_physical(omegaYFourier);
omegaZ = ft.transform_to_physical(omegaZFourier);
dzOmegaY = ft.transform_to_physical(dzOmegaYFourier);

% transform to plus units (velocity and vorticity transform in the same way)
uPlus = convert_to_plus_units(u, reBulk, reTau);
vPlus = convert_to_plus_units(v, reBulk, reTau);
wPlus = convert_to_plus_units(w, reBulk, reTau);
omegaXPlus = convert_to_plus_units(omegaX, reBulk, reTau);
omegaYPlus = convert_to_plus_units(omegaY, reBulk, reTau);
omegaZPlus = convert_to_plus_units(omegaZ, reBulk, reTau);
dzOmegaYPlus = convert_to_plus_units(dzOmegaY, reBulk, reTau);

% compute total stress
uMeanPlus = squeeze(mean(uPlus, [1 2]));
viscousStressPlus = dy.compute_derivative(uMeanPlus) / reTau;
uvReynoldsStressPlus = squeeze(mean(uPlus .* vPlus, [1 2]));
totalStressPlus = viscousStressPlus - uvReynoldsStressPlus;
slopeTotalStress = dy.compute_derivative(totalStressPlus);

% compute mean vorticity flux
vorticityFluxYZPlus = (vPlus .* omegaZPlus) - (wPlus .* omegaYPlus) ...
    - ((dy.compute_derivative(omegaZPlus) - dzOmegaYPlus) / reTau);
vorticityFluxYZMeanPlus = squeeze(mean(vorticityFluxYZPlus, [1 2]));

% plot quantities
figure;
plot(yCheb, totalStressPlus, 'k-', 'LineWidth', 1.5, 'DisplayName', 'total stress');
hold on;
plot(yCheb, viscousStressPlus, 'b-', 'LineWidth', 1.5, 'DisplayName', 'viscous stress');
plot(yCheb, -uvReynoldsStressPlus, 'r-', 'LineWidth', 1.5, 'DisplayName', '-1.0 x Reynolds stress');
xlabel('$y/h$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\tau^+$', 'Interpreter', 'latex', 'FontSize', 14)
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);

figure;
plot(yCheb, slopeTotalStress, 'k-', 'LineWidth', 1.5, 'DisplayName', '$\mathrm{d}_y \langle \tau_{tot}^+ \rangle$');
hold on;
plot(yCheb, vorticityFluxYZMeanPlus, 'r--', 'LineWidth', 1.5, 'DisplayName', '$\langle \Sigma_{yz}^+ \rangle$');
xlabel('$y/h$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$q h / u_{\tau}^2$', 'Interpreter', 'latex', 'FontSize', 14)
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 14);


function qPlus = convert_to_plus_units(qBulk, reBulk, reTau)
    qPlus = reBulk / reTau * qBulk;
end

