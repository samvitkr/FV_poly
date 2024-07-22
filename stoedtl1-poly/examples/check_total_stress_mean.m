% © 2023 Simon Toedtli <s.toedtli@jhu.edu>, all rights reserved
clear;  clc;  close all;
addpath '/home/skumar67/data-geyink1/skumar67/dmsuite'
% automatically add code modules to matlab path
functionCallStack = dbstack;
[scriptFolder, ~, ~] = fileparts(which(functionCallStack(1).file));
[sourceFolder, ~, ~] = fileparts(scriptFolder);
addpath(fullfile(sourceFolder, 'data_readers'), fullfile(sourceFolder, 'utilities'));

% user input: specify run parameters. This assumes your data is stored in a directory with full path dataFolder/runName
runName = 'newtonian_sheng_1';
fileNrStart = 49000;
fileNrEnd = 72000;
fileNrIncr = 1000;
ngx = 193;  % ng{i}: grid size of restart files in coordinate i, including periodic points in x, z.
ngy = 225;  %        ghost points do not have to be included here, they are accounted for in the data reading routines
ngz = 161;
nCheb = 175;  % number of Chebyshev collocation points in y (has to be determined a priori)
Lz = 2.0 * pi;
reBulk = 4667.0;
reTau = 253.0;
%dataFolder = fullfile(getenv('HOME'), 'research_data', 'polymer');

% read grid and generate operators
%runFolder = fullfile(dataFolder, runName);
runFolder='/home/skumar67/data-geyink1/skumar67/FV_wi_1/test_trans_wi1-0_1';
[xGridPointsDns, yGridPointsDns, ~, ~] = read_grid(runFolder, ngx, ngy);
deltaX = xGridPointsDns(2) - xGridPointsDns(1);
deltaZ = Lz / (ngz - 1);
zGridPointsDns = 0 : deltaZ : Lz;
[yCheb, ~] = chebdif(nCheb, 1);
interpY = InterpolationOperatorY(yGridPointsDns, yCheb);
ft = SpatialFourierTransform();

% pre-allocate variables
%uMean = zeros(ngy + 1, 1);
uvReynoldsStress = zeros(nCheb, 1);
uMean=uvReynoldsStress;
nSamples = 0;
m=matfile('mean_wi1.mat','Writable',true)
m. yGridPointsDns= yGridPointsDns;
for fileNr = fileNrStart : fileNrIncr : fileNrEnd
    % read DNS data: cartesian components in phyiscal domain, represented on staggered grid
    fileNrString = num2str(fileNr, '%07d');
    disp(fileNrString);
    [u, v, w] = read_velocity(runFolder, fileNrString, ngx, ngy, ngz);  % these arrays still contain ghost points
    [u, v, w] = remove_velocity_ghost_points(u, v, w);
    
    % compute mean velocity
    %uMean = uMean + squeeze(mean(u, [1 2]));

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

    % map back to physical space
    u = ft.transform_to_physical(uFourier);
    v = ft.transform_to_physical(vFourier);

	uMean = uMean + squeeze(mean(u, [1 2]));

    % compute Reynolds stress
    uvReynoldsStress = uvReynoldsStress + squeeze(mean(u .* v, [1 2]));
    nSamples = nSamples + 1;
end

% average in time
uMean = uMean / nSamples;
uvReynoldsStress = uvReynoldsStress / nSamples;

% compute to plus units
uMeanPlus = uMean * reBulk / reTau;
uvReynoldsStressPlus = uvReynoldsStress * (reBulk / reTau)^2;

% compute stresses
dy = FirstDerivativeYChebyshev(nCheb);
uMeanPlus = interpY.interpolate_u_to_chebyshev_grid(uMeanPlus);
viscousStressPlus = dy.compute_derivative(uMeanPlus) / reTau;
totalStressPlus = viscousStressPlus - uvReynoldsStressPlus;

% compute slopes
slopeViscousStress = dy.compute_derivative(viscousStressPlus);
slopeReynoldsStress = dy.compute_derivative(uvReynoldsStressPlus);
slopeTotalStress = dy.compute_derivative(totalStressPlus);


%m=matfile('mean_wi1.mat','Writable',true)
m.Um=uMean;
m.uvm=uvReynoldsStress;
m.yCheb=yCheb;

%% plot quantities
%yPlus = reTau * (yCheb + 1);
%figure;
%plot(yPlus, totalStressPlus, 'k-', 'LineWidth', 1.5, 'DisplayName', 'total stress');
%hold on;
%plot(yPlus, viscousStressPlus, 'b-', 'LineWidth', 1.5, 'DisplayName', 'viscous stress');
%plot(yPlus, -uvReynoldsStressPlus, 'r-', 'LineWidth', 1.5, 'DisplayName', '-1.0 x Reynolds stress');
%xlabel('$y^+$', 'Interpreter', 'latex', 'FontSize', 14);
%ylabel('$\tau^+$', 'Interpreter', 'latex', 'FontSize', 14)
%xlim([yPlus(end) yPlus(1)]);
%legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);
%title('time-averaged results');
%
%figure;
%plot(yPlus, slopeTotalStress, 'k-', 'LineWidth', 1.5);
%xlabel('$y^+$', 'Interpreter', 'latex', 'FontSize', 14);
%ylabel('$\mathrm{d}_y \langle \tau_{tot}^+ \rangle$', 'Interpreter', 'latex', 'FontSize', 14);
%xlim([yPlus(end) yPlus(1)]);
%title('time-averaged results');
%
%figure;
%plot(yPlus, slopeTotalStress, 'k-', 'LineWidth', 1.5, 'DisplayName', 'slope total stress');
%hold on;
%plot(yPlus, slopeViscousStress, 'b-', 'LineWidth', 1.5, 'DisplayName', 'slope viscous stress');
%plot(yPlus, - slopeReynoldsStress, 'r-', 'LineWidth', 1.5, 'DisplayName', '-1.0 x slope Reynolds stress');
%xlabel('$y^+$', 'Interpreter', 'latex', 'FontSize', 14);
%ylabel('$\mathrm{d}_y \langle \tau^+ \rangle$', 'Interpreter', 'latex', 'FontSize', 14);
%legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);
%xlim([yPlus(end) yPlus(1)]);
%title('time-averaged results');
%
%figure;
%semilogx(yPlus(2:end), slopeTotalStress(2:end), 'k-', 'LineWidth', 1.5);
%xlabel('$y^+$', 'Interpreter', 'latex', 'FontSize', 14);
%ylabel('$\mathrm{d}_y \langle \tau_{tot}^+ \rangle$', 'Interpreter', 'latex', 'FontSize', 14);
%xlim([1 reTau]);
%title('time-averaged results');
%
