% © 2023 Simon Toedtli <s.toedtli@jhu.edu>, all rights reserved
clear;  clc;  close all;
addpath '/home/skumar67/data-geyink1/skumar67/dmsuite'
sourceFolder='../';
% automatically add code modules to matlab path
functionCallStack = dbstack;
[scriptFolder, ~, ~] = fileparts(which(functionCallStack(1).file));
[sourceFolder, ~, ~] = fileparts(scriptFolder);
addpath(fullfile(sourceFolder, 'data_readers'), fullfile(sourceFolder, 'utilities'));

% user input: specify run parameters. This assumes your data is stored in a directory with full path dataFolder/runName
%runName = 'test_wi_6-67_1';
%fileNr = 79000;
ngx = 513;  % ng{i}: grid size of restart files in coordinate i, including periodic points in x, z.
ngy = 321;  %        ghost points do not have to be included here, they are accounted for in the data reading routines
ngz = 385;
nCheb = 220;  % number of Chebyshev collocation points in y (has to be determined a priori)
Lz = 2.0 * pi;
%dataFolder = fullfile(getenv('HOME'), 'research_data', 'polymer');
dataFolder ='/vast/geyink1/skumar67/FV_data/wi13p5'
Wi=13.5;
beta=0.9;
L_max=100;
a=1-3/L_max^2;
re=4667
% read grid and generate operators
%runFolder = fullfile(dataFolder, runName);
%runFolder= 'C:\Users\samvi\Dropbox\SimonsProject\Finite_Vol_data\Visco_data_sample';
%runFolder='/home/skumar67/data-geyink1/skumar67/FV_visco';
runFolder='/home/skumar67/data-geyink1/skumar67/FV_poly_codes/FV_wi_13p5/data'
%runFolder='/home/skumar67/data-geyink1/skumar67/FV_data_'

[xGridPointsDns, yGridPointsDns, ~, ~] = read_grid(runFolder, ngx, ngy);
deltaX = xGridPointsDns(2) - xGridPointsDns(1);
deltaZ = Lz / (ngz - 1);
[yCheb, ~] = chebdif(nCheb, 1);
interpY = InterpolationOperatorY(yGridPointsDns, yCheb);
ft = SpatialFourierTransform();

%piox=zeros(ngz-1,ngx-1,nCheb);
%pioy=zeros(ngz-1,ngx-1,nCheb);
%pioz=zeros(ngz-1,ngx-1,nCheb);
%pioxF=zeros(ngz-1,ngx-1,nCheb);
%pioyF=zeros(ngz-1,ngx-1,nCheb);
%piozF=zeros(ngz-1,ngx-1,nCheb);


fstart=0000;
fend=240000;
fskip=1000;
nf=(fend-fstart)/fskip+1;
for fileNr=fstart:fskip:fend
fileNr

% read DNS data: Cartesian velocities and confirmation tensor in physical domain, represented on staggered grid
fileNrString = num2str(fileNr, '%07d');
%[u, v, w] = read_velocity(runFolder, fileNrString, ngx, ngy, ngz);  % these arrays still contain ghost points
confTensor = read_confirmation_tensor(dataFolder, fileNrString, ngx, ngy, ngz);
%[u, v, w] = remove_velocity_ghost_points(u, v, w);
confTensor = remove_confirmation_tensor_ghost_points(confTensor);

[confTensorFourier,kx,kz] = ft.transform_confirmation_tensor_to_fourier(confTensor,deltaX, deltaZ);
size(kz)
size(kz)
% collocate data spectrally in x, z
shiftX = -0.5 * deltaX;
shiftZ = -0.5 * deltaZ;
confTensorFourier = interpolate_confirmation_tensor_fourier_in_x(confTensorFourier, kx, shiftX);
confTensorFourier = interpolate_confirmation_tensor_fourier_in_z(confTensorFourier, kz, shiftZ);


dx = FirstDerivativeXFourier(kx.get_radial_frequency_vector());
dz = FirstDerivativeZFourier(kz.get_radial_frequency_vector());
dy = FirstDerivativeYChebyshev(nCheb);


confTensorFourier = interpY.interpolate_confirmation_tensor_to_chebyshev_grid(confTensorFourier);
confTensorPhysical=ft.transform_confirmation_tensor_to_physical(confTensorFourier);
psi=1-(confTensorPhysical.Cxx+confTensorPhysical.Cyy+confTensorPhysical.Czz)./L_max^2;

polystress.tauxx=(confTensorPhysical.Cxx./psi - 1/a)./Wi;
polystress.tauyy=(confTensorPhysical.Cyy./psi - 1/a)./Wi;
polystress.tauzz=(confTensorPhysical.Czz./psi - 1/a)./Wi;
polystress.tauxy=(confTensorPhysical.Cxy./psi )./Wi;
polystress.tauxz=(confTensorPhysical.Cxz./psi )./Wi;
polystress.tauyz=(confTensorPhysical.Cyz./psi )./Wi;

polystressFourier.tauxx=ft.transform_to_fourier(polystress.tauxx);
polystressFourier.tauxz=ft.transform_to_fourier(polystress.tauxz);
polystressFourier.tauxy=ft.transform_to_fourier(polystress.tauxy);
polystressFourier.tauyz=ft.transform_to_fourier(polystress.tauyz);
polystressFourier.tauzz=ft.transform_to_fourier(polystress.tauzz);

dxfxx=ft.transform_to_physical(dx.compute_derivative(polystressFourier.tauxx));
dyfxy=dy.compute_derivative(polystress.tauxy);
dzfxz=ft.transform_to_physical(dz.compute_derivative(polystressFourier.tauxz));

dxfxy=ft.transform_to_physical(dx.compute_derivative(polystressFourier.tauxy));
dyfyy=dy.compute_derivative(polystress.tauyy);
dzfyz=ft.transform_to_physical(dz.compute_derivative(polystressFourier.tauyz));

dxfxz=ft.transform_to_physical(dx.compute_derivative(polystressFourier.tauxz));
dyfyz=dy.compute_derivative(polystress.tauyz);
dzfzz=ft.transform_to_physical(dz.compute_derivative(polystressFourier.tauzz));

polyx=((1-beta)./re)*(dxfxx+dyfxy+dzfxz);
polyy=((1-beta)./re)*(dxfxy+dyfyy+dzfyz);
polyz=((1-beta)./re)*(dxfxz+dyfyz+dzfzz);

fnt=sprintf("transferfields_%07d.mat",fileNr)
fmt=fullfile(runFolder,fnt)
mt=matfile(fmt,'Writable',true)
mt.polyx=polyx;
mt.polyy=polyy;
mt.polyz=polyz;

%
%
%polyxF=ft.transform_to_fourier(polyx);
%polyyF=ft.transform_to_fourier(polyy);
%polyzF=ft.transform_to_fourier(polyz);


%torquex=			 dy.compute_derivative(polyz)  - ft.transform_to_physical(dz.compute_derivative(polyyF));
%torquey=ft.transform_to_physical(dz.compute_derivative(polyxF))- ft.transform_to_physical(dx.compute_derivative(polyzF));
%torquez=ft.transform_to_physical(dx.compute_derivative(polyyF))-			  dy.compute_derivative(polyx);
%
%torquexF=ft.transform_to_fourier(torquex);
%torqueyF=ft.transform_to_fourier(torquey);
%torquezF=ft.transform_to_fourier(torquez);

%fng=sprintf("velgrad_%07d.mat",fileNr)
%fmg=fullfile(runFolder,fng);
%mg=matfile(fmg);
%omegaxF=ft.transform_to_fourier(mg.dwdy-mg.dvdz);
%omegayF=ft.transform_to_fourier(mg.dudz-mg.dwdx);
%omegazF=ft.transform_to_fourier(mg.dvdx-mg.dudy);
%
%omegax=(mg.dwdy-mg.dvdz);
%omegay=(mg.dudz-mg.dwdx);
%omegaz=(mg.dvdx-mg.dudy);

%pioxF=pioxF+omegaxF.*conj( torquexF );
%pioyF=pioyF+omegayF.*conj( torqueyF );
%piozF=piozF+omegazF.*conj( torquezF );

%piox=piox+omegax.*torquex;
%pioy=pioy+omegay.*torquey;
%pioz=pioz+omegaz.*torquez;
%
%
end


%fme=fullfile(runFolder,'enstrophy_transferF.mat')
%fme=fullfile(runFolder,'enstrophy_transfer.mat')
%me=matfile(fme,'Writable',true)
%me.piox=piox./nf;
%me.pioy=pioy./nf;
%me.pioz=pioz./nf;
%me.pioxF=pioxF./nf;
%me.pioyF=pioyF./nf;
%me.piozF=piozF./nf;
%fny="ygrid.mat";
% % %fmy="/home/skumar67/data-geyink1/skumar67/FV_visco/ygrid.mat";
%fmy=fullfile(runFolder,fny)
%my=matfile(fmy,'Writable',true)
%my.yCheb=yCheb;



