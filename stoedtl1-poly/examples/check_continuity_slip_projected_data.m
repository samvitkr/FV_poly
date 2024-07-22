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
[xGridPointsDns, yGridPointsDns, ~, yCellCenter] = read_grid(runFolder, ngx, ngy);
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

% save fields that are collocated in x, z
uCollocatedXZ = ft.transform_to_physical(uFourier);
wCollocatedXZ = ft.transform_to_physical(wFourier);

% collocate on Chebyshev grid in y using spline interpolation
uFourier = interpY.interpolate_u_to_chebyshev_grid(uFourier);
vFourier = interpY.interpolate_v_to_chebyshev_grid(vFourier);
wFourier = interpY.interpolate_w_to_chebyshev_grid(wFourier);

% the velocity fields were divergence-free w.r.t. the finite volume numerics, but not necessarily w.r.t the spectral
% numerics. We therefore project onto a divergence free subspace (note: the projected fields satisfy no-through at the
% wall, but they may have non-zero slip).
[uFourier, vFourier, wFourier] = make_field_divergence_free(uFourier, vFourier, wFourier, kx, kz, nCheb);

% do a number of sanity checks on the data
% 1) find the maximum continuity residual as a function of y
dx = FirstDerivativeXFourier(kx.get_radial_frequency_vector());
dz = FirstDerivativeZFourier(kz.get_radial_frequency_vector());
dy = FirstDerivativeYChebyshev(nCheb);
residualFourier = dx.compute_derivative(uFourier) + dy.compute_derivative(vFourier) + dz.compute_derivative(wFourier);
residual = ft.transform_to_physical(residualFourier);
maxResidualAtEachY = squeeze(max(abs(residual), [], [1 2]));

rgbOrange = [0.850, 0.325, 0.098];
figure();
subplot(2, 1, 1);
plot(yCheb, maxResidualAtEachY, 'Color', rgbOrange, 'LineWidth', 2.0);
hold on;
text(-0.45, 0.35, 'Including wall points', 'FontSize', 14);
xlabel('$y/h$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\| \partial_x u + \partial_y v + \partial_z w  \|_{\infty}$', 'Interpreter', 'latex', 'FontSize', 14);
subplot(2, 1, 2);
plot(yCheb(2:end-1), maxResidualAtEachY(2:end-1), 'Color', rgbOrange, 'LineWidth', 2.0);
hold on;
text(-0.45, 1.25e-12, 'Excluding wall points', 'FontSize', 14);
xlabel('$y/h$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\| \partial_x u + \partial_y v + \partial_z w  \|_{\infty}$', 'Interpreter', 'latex', 'FontSize', 14);

% 2) plot slip velocity at bottom wall
uProjected = ft.transform_to_physical(uFourier);
vProjected = ft.transform_to_physical(vFourier);
wProjected = ft.transform_to_physical(wFourier);
idxYBot = nCheb;  % yCheb goes from +1 to -1

figure();  % u
contourf(xGridPointsDns(1:end-1), zGridPointsDns(1:end-1), uProjected(:, :, idxYBot));
colorbar();
xlabel('$x/h$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$z/h$', 'Interpreter', 'latex', 'FontSize', 14);
title('Interpolated and projected: u(y_w) / U_b');

figure();  % v
contourf(xGridPointsDns(1:end-1), zGridPointsDns(1:end-1), vProjected(:, :, idxYBot));
colorbar();
xlabel('$x/h$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$z/h$', 'Interpreter', 'latex', 'FontSize', 14);
title('Interpolated and projected: v(y_w) / U_b');

figure();  % w
contourf(xGridPointsDns(1:end-1), zGridPointsDns(1:end-1), wProjected(:, :, idxYBot));
colorbar();
xlabel('$x/h$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$z/h$', 'Interpreter', 'latex', 'FontSize', 14);
title('Interpolated and projected: w(y_w) / U_b');

% 3) compare the wall-normal derivative of u and w at the wall with a reference value
% derivatives from projected data
dudyProjected = dy.compute_derivative(uProjected);
dwdyProjected = dy.compute_derivative(wProjected);
dudyProjectedBot = dudyProjected(:, :, idxYBot);
dwdyProjectedBot = dwdyProjected(:, :, idxYBot);
meanWallShearStressProjected = mean(dudyProjectedBot, 'all');
% reference data
idxYBotFiniteDiff = 2;  % finite difference grid is 0 -> 2. The first index corresponds to the ghost point
deltaY = yCellCenter(1);
deltaU = uCollocatedXZ(:, :, idxYBotFiniteDiff);
deltaW = wCollocatedXZ(:, :, idxYBotFiniteDiff);
dudyFiniteDiffBot = deltaU / deltaY;
dwdyFiniteDiffBot = deltaW / deltaY;
meanWallShearStressFiniteDiff = mean(dudyFiniteDiffBot, 'all');
% compute normalized difference: in both cases the plots are normalized w.r.t. the scale that makes the finite
% difference wall shear stress spectrum O(1)
normalizedDifferenceDudy = abs(dudyProjectedBot - dudyFiniteDiffBot) / max(abs(dudyFiniteDiffBot), [], 'all');
normalizedDifferenceDwdy = abs(dwdyProjectedBot - dwdyFiniteDiffBot) / max(abs(dwdyFiniteDiffBot), [], 'all');
differenceMeanWallShearStress = abs(meanWallShearStressProjected - meanWallShearStressFiniteDiff);
% compare the mean wall-shear stress
disp(['instantaneous plane averaged wall shear stress spectral: ' num2str(meanWallShearStressProjected)]);
disp(['instantaneous plane averaged wall shear stress finite diff: ' num2str(meanWallShearStressFiniteDiff)]);
disp(['absolute difference = ' num2str(differenceMeanWallShearStress)]);
disp(['relative difference = ' num2str(abs(differenceMeanWallShearStress / meanWallShearStressFiniteDiff))]);

figure();  % dudy
contourf(xGridPointsDns(1:end-1), zGridPointsDns(1:end-1), normalizedDifferenceDudy);
colorbar();
xlabel('$x/h$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$z/h$', 'Interpreter', 'latex', 'FontSize', 14);
title('Normalized difference FD vs spectral: du/dy(y_w)');

figure();  % dwdy
contourf(xGridPointsDns(1:end-1), zGridPointsDns(1:end-1), normalizedDifferenceDwdy);
colorbar();
xlabel('$x/h$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$z/h$', 'Interpreter', 'latex', 'FontSize', 14);
title('Normalized difference finite difference and spectral: dw/dy(y_w)');
