% © 2022 Simon Toedtli <s.toedtli@jhu.edu>, all rights reserved
clear;  clc;  close all;

% automatically add code modules to matlab path
functionCallStack = dbstack;
[scriptFolder, ~, ~] = fileparts(which(functionCallStack(1).file));
[sourceFolder, ~, ~] = fileparts(scriptFolder);
addpath(fullfile(sourceFolder, 'data_readers'), fullfile(sourceFolder, 'utilities'));

% user input: specify run parameters. This assumes your data is stored in a directory with full path dataFolder/runName
runName = 'canonical_4pi_2pi_5';
ngx = 193;  % ng{i}: grid size of restart files in coordinate i, including periodic points in x, z.
ngy = 225;  %        ghost points do not have to be included here, they are accounted for in the data reading routines
ngz = 161;
nChebVector = 10 : 5 : ngy;  % vector containing the number of collocation points you want to test
fileNr = 400;
idxX = 34;   % idx{X,Z} specify a sample x, z location in the data where you want to test the reconstruction
idxZ = 120;  % to get a more robust estimate, you may want to average over several (x, z) locations
dataFolder = '/home/simon/research_data';

% load data
runFolder = fullfile(dataFolder, runName);
fileNrString = num2str(fileNr, '%07d');
[~, yGridPointsDns, ~, ~] = read_grid(runFolder, ngx, ngy);
[u, v, w] = read_velocity(runFolder, fileNrString, ngx, ngy, ngz);
[u, v, w] = remove_velocity_ghost_points(u, v, w);
vTruth = squeeze(v(idxZ, idxX, :));

% construct Chebyshev grid and interpolate
nTestResolutions = numel(nChebVector);
absoluteErrorTwoNorm = zeros(1, nTestResolutions);
relativeErrorTwoNorm = zeros(1, nTestResolutions);
normTruth = norm(vTruth, 2);
for i = 1 : nTestResolutions
    nCheb = nChebVector(i);
    disp(nCheb);
    [yCheb, ~] = chebdif(nCheb, 1);
    interpY = InterpolationOperatorY(yGridPointsDns, yCheb);
    vCollocation = interpY.interpolate_v_to_chebyshev_grid(vTruth);  % DNS grid -> collocation points using spline interpolation
    % collocation points -> DNS grid using Chebyshev interpolant. We have to use interpY.get_y_grid_points_dns() instead of
    % yGridPointsDns, because only the former is over the interval [-1, 1]
    vReconstr = chebint(vCollocation, interpY.get_y_grid_points_dns());  % collocation points -> DNS grid using Chebyshev interpolant
    absoluteError = abs(vReconstr - vTruth);
    absoluteErrorTwoNorm(i) = norm(absoluteError, 2);
    relativeErrorTwoNorm(i) = norm(absoluteError, 2) / normTruth;
end

fsLabel = 16;
figure();  % plot absolute error
semilogy(nChebVector, absoluteErrorTwoNorm, 'ks')
xlabel('$N$', 'Interpreter', 'latex', 'FontSize', fsLabel);
ylabel('$\| v - P_N(S_N(v)) \|_2$', 'Interpreter', 'latex', 'FontSize', fsLabel);
title('Absolute reconstruction error');

figure();  % plot relative error
semilogy(nChebVector, relativeErrorTwoNorm, 'ks')
xlabel('$N$', 'Interpreter', 'latex', 'FontSize', fsLabel);
ylabel('$\| v - P_N(S_N(v)) \|_2 / \|v \|_2$', 'Interpreter', 'latex', 'FontSize', fsLabel);
title('Relative reconstruction error');

rgbGray = 0.8 * [1, 1, 1];
figure();  % plot example reconstructed wall-normal profile
plot(interpY.get_y_grid_points_dns(), vTruth, 'LineStyle', '-', 'Color', rgbGray, 'LineWidth', 2.0, 'DisplayName', 'original');
hold on;
plot(interpY.get_y_grid_points_dns(), vReconstr, 'k-', 'LineWidth', 2.0, 'DisplayName', 'reconstructed');
xlabel('y/h');
ylabel('v(x_0, y, z_0)');
legend('Location', 'best');
