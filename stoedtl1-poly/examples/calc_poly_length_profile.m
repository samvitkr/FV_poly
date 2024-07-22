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
dataFolder = fullfile(getenv('HOME'), 'research_data', 'polymer');
Wi=13.5;
beta=0.9;
L_max=100;
a=1-3/L_max^2;
re=4667
% read grid and generate operators
%runFolder = fullfile(dataFolder, runName);
%runFolder= 'C:\Users\samvi\Dropbox\SimonsProject\Finite_Vol_data\Visco_data_sample';
%runFolder='/home/skumar67/data-geyink1/skumar67/FV_visco';
runFolder='/home/skumar67/data-geyink1/skumar67/FV_wi_13p5'
%runFolder='/home/skumar67/data-geyink1/skumar67/FV_data_'

[xGridPointsDns, yGridPointsDns, ~, ~] = read_grid(runFolder, ngx, ngy);
deltaX = xGridPointsDns(2) - xGridPointsDns(1);
deltaZ = Lz / (ngz - 1);
[yCheb, ~] = chebdif(nCheb, 1);
interpY = InterpolationOperatorY(yGridPointsDns, yCheb);
ft = SpatialFourierTransform();


fstart=0000;
fend=240000;
fskip=1000;

for fileNr=fstart:fskip:fend
fileNr

% read DNS data: Cartesian velocities and confirmation tensor in physical domain, represented on staggered grid
fileNrString = num2str(fileNr, '%07d');
%[u, v, w] = read_velocity(runFolder, fileNrString, ngx, ngy, ngz);  % these arrays still contain ghost points
confTensor = read_confirmation_tensor(runFolder, fileNrString, ngx, ngy, ngz);
%[u, v, w] = remove_velocity_ghost_points(u, v, w);
confTensor = remove_confirmation_tensor_ghost_points(confTensor);

% tranform to Fourier domain
%[uFourier, kx, kz] = ft.transform_to_fourier(u, deltaX, deltaZ);
%vFourier = ft.transform_to_fourier(v);  % no need to pass sampling rate, wavenumber object is needed only once
%wFourier = ft.transform_to_fourier(w);
[confTensorFourier,kx,kz] = ft.transform_confirmation_tensor_to_fourier(confTensor,deltaX, deltaZ);
size(kz)
size(kz)
% collocate data spectrally in x, z
shiftX = -0.5 * deltaX;
shiftZ = -0.5 * deltaZ;
%uFourier = interpolate_fourier_in_z(uFourier, kz, shiftZ);
%vFourier = interpolate_fourier_in_x(vFourier, kx, shiftX);
%vFourier = interpolate_fourier_in_z(vFourier, kz, shiftZ);
%wFourier = interpolate_fourier_in_x(wFourier, kx, shiftX);
confTensorFourier = interpolate_confirmation_tensor_fourier_in_x(confTensorFourier, kx, shiftX);
confTensorFourier = interpolate_confirmation_tensor_fourier_in_z(confTensorFourier, kz, shiftZ);

% collocate on Chebyshev grid in y using spline interpolation
%uFourier = interpY.interpolate_u_to_chebyshev_grid(uFourier);
%vFourier = interpY.interpolate_v_to_chebyshev_grid(vFourier);
%wFourier = interpY.interpolate_w_to_chebyshev_grid(wFourier);

dx = FirstDerivativeXFourier(kx.get_radial_frequency_vector());
dz = FirstDerivativeZFourier(kz.get_radial_frequency_vector());
dy = FirstDerivativeYChebyshev(nCheb);


confTensorFourier = interpY.interpolate_confirmation_tensor_to_chebyshev_grid(confTensorFourier);
confTensorPhysical=ft.transform_confirmation_tensor_to_physical(confTensorFourier);
stretch=(confTensorPhysical.Cxx+confTensorPhysical.Cyy+confTensorPhysical.Czz)./L_max^2;
cxxl=confTensorPhysical.Cxx./L_max^2;
cyyl=confTensorPhysical.Cyy./L_max^2;
czzl=confTensorPhysical.Czz./L_max^2;
cxyl=confTensorPhysical.Cxy./L_max^2;
%polystress.tauxx=(confTensorPhysical.Cxx./psi - 1/a)./Wi;
%polystress.tauyy=(confTensorPhysical.Cyy./psi - 1/a)./Wi;
%polystress.tauzz=(confTensorPhysical.Czz./psi - 1/a)./Wi;
%polystress.tauxy=(confTensorPhysical.Cxy./psi )./Wi;
%polystress.tauxz=(confTensorPhysical.Cxz./psi )./Wi;
%polystress.tauyz=(confTensorPhysical.Cyz./psi )./Wi;
%
%polystressFourier.tauxx=ft.transform_to_fourier(polystress.tauxx);
%polystressFourier.tauxz=ft.transform_to_fourier(polystress.tauxz);
%polystressFourier.tauxy=ft.transform_to_fourier(polystress.tauxy);
%polystressFourier.tauyz=ft.transform_to_fourier(polystress.tauyz);
%polystressFourier.tauzz=ft.transform_to_fourier(polystress.tauzz);


%fxxF=dx.compute_derivative(polystressFourier.tauxx);
%fxzF=dz.compute_derivative(polystressFourier.tauxz);

%dxfxx=ft.transform_to_physical(dx.compute_derivative(polystressFourier.tauxx));
%dyfxy=dy.compute_derivative(polystress.tauxy);
%dzfxz=ft.transform_to_physical(dz.compute_derivative(polystressFourier.tauxz));


%dxfxy=ft.transform_to_physical(dx.compute_derivative(polystressFourier.tauxy));
%dyfyy=dy.compute_derivative(polystress.tauyy);
%dzfyz=ft.transform_to_physical(dz.compute_derivative(polystressFourier.tauyz));
%
%
%dxfxz=ft.transform_to_physical(dx.compute_derivative(polystressFourier.tauxz));
%dyfyz=dy.compute_derivative(polystress.tauyz);
%dzfzz=ft.transform_to_physical(dz.compute_derivative(polystressFourier.tauzz));

%polyx=((1-beta)./re)*(dxfxx+dyfxy+dzfxz);
%polyy=((1-beta)./re)*(dxfxy+dyfyy+dzfyz);
%polyz=((1-beta)./re)*(dxfxz+dyfyz+dzfzz);
%
%
%polyxF=ft.transform_to_fourier(polyx);
%polyyF=ft.transform_to_fourier(polyy);
%polyzF=ft.transform_to_fourier(polyz);
%
%torque_x=dy.compute_derivative(polyz)-ft.transform_to_physical(dz.compute_derivative( polyyF ));
%torque_y=ft.transform_to_physical(dz.compute_derivative( polyxF ))-ft.transform_to_physical(dx.compute_derivative( polyzF ));
%torque_z=ft.transform_to_physical(dx.compute_derivative( polyyF ))-dy.compute_derivative(polyx);
%

% the velocity fields were divergence-free w.r.t. the finite volume numerics, but not necessarily w.r.t the spectral
% numerics. We therefore project onto a divergence free subspace (note: the projected fields satisfy no-through at the
% wall, but they may have non-zero slip).
%[uFourier, vFourier, wFourier] = make_field_divergence_free(uFourier, vFourier, wFourier, kx, kz, nCheb);
%
%dudxF=(dx.compute_derivative(uFourier));
%dvdxF=(dx.compute_derivative(vFourier));
%dwdxF=(dx.compute_derivative(wFourier));
%dudyF=(dy.compute_derivative(uFourier));
%dvdyF=(dy.compute_derivative(vFourier));
%dwdyF=(dy.compute_derivative(wFourier));
%dudzF=(dz.compute_derivative(uFourier));
%dvdzF=(dz.compute_derivative(vFourier));
%dwdzF=(dz.compute_derivative(wFourier));
%
%dudx=ft.transform_to_physical(dx.compute_derivative(uFourier));
%dvdx=ft.transform_to_physical(dx.compute_derivative(vFourier));
%dwdx=ft.transform_to_physical(dx.compute_derivative(wFourier));
%dudy=ft.transform_to_physical(dy.compute_derivative(uFourier));
%dvdy=ft.transform_to_physical(dy.compute_derivative(vFourier));
%dwdy=ft.transform_to_physical(dy.compute_derivative(wFourier));
%dudz=ft.transform_to_physical(dz.compute_derivative(uFourier));
%dvdz=ft.transform_to_physical(dz.compute_derivative(vFourier));
%dwdz=ft.transform_to_physical(dz.compute_derivative(wFourier));
%
%
%
%omegaXFourier = dy.compute_derivative(wFourier) - dz.compute_derivative(vFourier);
%omegaYFourier = dz.compute_derivative(uFourier) - dx.compute_derivative(wFourier);
%omegaZFourier = dx.compute_derivative(vFourier) - dy.compute_derivative(uFourier);
%omegaX = ft.transform_to_physical(omegaXFourier);
%omegaY = ft.transform_to_physical(omegaYFourier);
%omegaZ = ft.transform_to_physical(omegaZFourier);
%
%viscF=  dx.compute_derivative( dx.compute_derivative(uFourier))+...
%        dy.compute_derivative( dy.compute_derivative(uFourier))+...
%        dz.compute_derivative( dz.compute_derivative(uFourier));
%visc=beta*ft.transform_to_physical(viscF)./re;
%
%ufield=ft.transform_to_physical(uFourier);
%vfield=ft.transform_to_physical(vFourier);
%wfield=ft.transform_to_physical(wFourier);
%
%voz=vfield.*omegaZ;
%woy=wfield.*omegaY;
%
%vozFourier=vFourier.*conj(omegaZFourier);
%woyFourier=wFourier.*conj(omegaYFourier);



% %fm=sprintf("/home/skumar67/data-geyink1/skumar67/FV_visco/velfields_%07d.mat",fileNr)
%fnm=sprintf("velfields_%07d.mat",fileNr);
%fm=fullfile(runFolder,fnm)
%m=matfile(fm,'Writable',true)
%m.uf=u;
%m.vf=v;
%m.wf=w;
%m.uFourier=uFourier;
%m.vFourier=vFourier;
%m.wFourier=wFourier;
%m.ufield=ft.transform_to_physical(uFourier);
%m.vfield=ft.transform_to_physical(vFourier);
%m.wfield=ft.transform_to_physical(wFourier);
%
%
%
% %fmo=sprintf("/home/skumar67/data-geyink1/skumar67/Visco_data_sample/vortfields_%07d.mat",fileNr)
% %fmo=sprintf("/home/skumar67/data-geyink1/skumar67/FV_visco/vortfields_%07d.mat",fileNr)
%fno=sprintf("vortfields_%07d.mat",fileNr);
%fmo=fullfile(runFolder,fno)
%mo=matfile(fmo,'Writable',true)
%mo.omegaXFourier=omegaXFourier;
%mo.omegaYFourier=omegaYFourier;
%mo.omegaZFourier=omegaXFourier;
%mo.omegaX=omegaX;
%mo.omegaY=omegaY;
%mo.omegaZ=omegaZ;
%
%fmt=sprintf("/home/skumar67/data-geyink1/skumar67/Visco_data_sample/transferfields_%07d.mat",fileNr)
fnt=sprintf("transferfields_%07d.mat",fileNr)
fmt=fullfile(runFolder,fnt)
mt=matfile(fmt,'Writable',true)
mt.l=stretch;
mt.cxxl=cxxl;
mt.cyyl=cyyl;
mt.czzl=czzl;
mt.cxyl=cxyl;
%mt.voz=voz;
%mt.woy=woy;
%mt.visc=visc;
%mt.vozF=vozFourier;
%mt.woyF=woyFourier;
%mt.poly=((1-beta)./re)*(fxx+fxy+fxz);
%mt.poly=polyx

%fmp=sprintf("/home/skumar67/data-geyink1/skumar67/Visco_data_sample/polyforce_%07d.mat",fileNr)
%fnp=sprintf("polyfields_%07d.mat",fileNr)
%fmp=fullfile(runFolder,fnp)
%mp=matfile(fmp,'Writable',true)
% %mp.dxfxx=fxx;
% %mp.dyfxy=fxy;
% %mp.dzfxz=fxz;
%mp.polyx=polyx;%((1-beta)./re)*(dxfxx+dyfxy+dzfxz);
%mp.polyy=polyy;%((1-beta)./re)*(dxfxy+dyfyy+dzfyz);
%mp.polyz=polyz;%((1-beta)./re)*(dxfxz+dyfyz+dzfzz);
%mp.torque_x=torque_x;
%mp.torque_y=torque_y;
%mp.torque_z=torque_z;

%fng=sprintf("velgrad_%07d.mat",fileNr)
%fmg=fullfile(runFolder,fng)
%mg=matfile(fmg,'Writable',true)
%
%mg.dudx=dudx;
%mg.dvdx=dvdx;
%mg.dwdx=dwdx;
%mg.dudy=dudy;
%mg.dvdy=dvdy;
%mg.dwdy=dwdy;
%mg.dudz=dudz;
%mg.dvdz=dvdz;
%mg.dwdz=dwdz;
%
%mg.dudxF=dudxF;
%mg.dvdxF=dvdxF;
%mg.dwdxF=dwdxF;
%mg.dudyF=dudyF;
%mg.dvdyF=dvdyF;
%mg.dwdyF=dwdyF;
%mg.dudzF=dudzF;
%mg.dvdzF=dvdzF;
%mg.dwdzF=dwdzF;
%
%mg.kx=kx;
%mg.kz=kz;
end




%fny="ygrid.mat";
% % %fmy="/home/skumar67/data-geyink1/skumar67/FV_visco/ygrid.mat";
%fmy=fullfile(runFolder,fny)
%my=matfile(fmy,'Writable',true)
%my.yCheb=yCheb;



