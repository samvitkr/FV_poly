clear;  clc;  close all;
addpath('../data_readers/');

runName = 'test_wi_6-67_1';
fileNr = 10;
ngx = 193;
ngy = 225;
ngz = 161;
nCheb = 151;
Lz = 2.0 * pi;
dataFolder = fullfile(getenv('HOME'), 'research_data', 'polymer');

% read grid and generate interpolation object
runFolder = fullfile(dataFolder, runName);
[~, yGridPointsDns, ~, ~] = read_grid(runFolder, ngx, ngy);
% deltaX = xGridPointsDns(2) - xGridPointsDns(1);
% deltaZ = Lz / (ngz - 1);
[yCheb, ~] = chebdif(nCheb, 1);
interpY = InterpolationOperatorY(yGridPointsDns, yCheb);

% read velocity field and generate a mock-up confirmation tensor to test the new interpolation routines
fileNrString = num2str(fileNr, '%07d');
[u, v, w] = read_velocity(runFolder, fileNrString, ngx, ngy, ngz);
[u, v, w] = remove_velocity_ghost_points(u, v, w);
confirmationTensor.Cxx = u;
confirmationTensor.Cyy = u;
confirmationTensor.Czz = u;
confirmationTensor.Cxy = u;
confirmationTensor.Cxz = u;
confirmationTensor.Cyz = u;

% interpolate in y and check error. The no-slip boundary condition, which is enforced in the u interpolation, should be
% recovered up to machine precision from the confirmation tensor interpolation
uCheb = interpY.interpolate_u_to_chebyshev_grid(u);
confirmationTensorCheb = interpY.interpolate_confirmation_tensor_to_chebyshev_grid(confirmationTensor);
diffXX = confirmationTensorCheb.Cxx - uCheb;
display_global_max(diffXX, 'Cxx');
diffYY = confirmationTensorCheb.Cyy - uCheb;
display_global_max(diffYY, 'Cyy');
diffZZ = confirmationTensorCheb.Czz - uCheb;
display_global_max(diffZZ, 'Czz');
diffXY = confirmationTensorCheb.Cxy - uCheb;
display_global_max(diffXY, 'Cxy');
diffXZ = confirmationTensorCheb.Cxz - uCheb;
display_global_max(diffXZ, 'Cxz');
diffYZ = confirmationTensorCheb.Cyz - uCheb;
display_global_max(diffYZ, 'Cyz');


function qMax = get_global_max(q)
    qMax = max(abs(q), [], 'all');
end

function display_global_max(q, qStr)
    disp(['Maximum deviation in ', qStr, ' = ', num2str(get_global_max(q))]);
end
