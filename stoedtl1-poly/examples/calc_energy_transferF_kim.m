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

ep=zeros(ngz-1,ngx-1,nCheb);
epF=zeros(ngz-1,ngx-1,nCheb);
pix=ep;
piy=ep;
piz=ep;
pixF=epF;
piyF=epF;
pizF=epF;


dis=zeros(ngz-1,ngx-1,nCheb);
fstart=0000;
fend=240000;
fskip=1000;
nf=(fend-fstart)/fskip+1;
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


%tau11F=ft.transform_to_fourier(polystress.tauxx);
%tau13F=ft.transform_to_fourier(polystress.tauxz);
%tau12F=ft.transform_to_fourier(polystress.tauxy);
%tau23F=ft.transform_to_fourier(polystress.tauyz);
%tau33F=ft.transform_to_fourier(polystress.tauzz);
%tau22F=ft.transform_to_fourier(polystress.tauyy);

fng=sprintf("velfields_%07d.mat",fileNr)
fmg=fullfile(runFolder,fng)
mg=matfile(fmg)

pixinst=polyx.*mg.ufield;
piyinst=polyy.*mg.vfield;
pizinst=polyz.*mg.wfield;

pixFinst=(ft.transform_to_fourier(polyx)).*conj(mg.uFourier);
piyFinst=(ft.transform_to_fourier(polyy)).*conj(mg.vFourier);
pizFinst=(ft.transform_to_fourier(polyz)).*conj(mg.wFourier);

%S11F=ft.transform_to_fourier(mg.dudx);
%S22F=ft.transform_to_fourier(mg.dvdy);
%S33F=ft.transform_to_fourier(mg.dwdz);
%S12F=ft.transform_to_fourier(0.5*( mg.dudy+mg.dvdx));
%S13F=ft.transform_to_fourier(0.5*( mg.dudz+mg.dwdx));
%S23F=ft.transform_to_fourier(0.5*( mg.dvdz+mg.dwdy));
%
%
%pixFinst=((1-beta)/re)*(tau11F.*conj(S11F)+...
%		        tau12F.*conj( ft.transform_to_fourier(mg.dudy)) + ...
%			tau13F.*conj( ft.transform_to_fourier(mg.dudz)));
%
%piyFinst=((1-beta)/re)*(tau12F.*conj(ft.transform_to_fourier(mg.dvdx))+...
%                        tau22F.*conj(S22F) +...
%			tau23F.*conj(ft.transform_to_fourier(mg.dvdz)));
%
%pizFinst=((1-beta)/re)*(tau13F.*conj(ft.transform_to_fourier(mg.dwdx))+...
%			tau23F.*conj(ft.transform_to_fourier(mg.dwdy)) +...
%		       	tau33F.*conj(S33F));
%
%
%pixinst=((1-beta)/re)* (polystress.tauxx.*(S11F)+...
%                        polystress.tauxy.*((mg.dudy)) + ...
%                        polystress.tauxz.*((mg.dudz)));
%
%piyinst=((1-beta)/re)* (polystress.tauxy.*((mg.dvdx))+...
%                        polystress.tauyy.*(S22F) +...
%                        polystress.tauyz.*((mg.dvdz)));
%
%pizinst=((1-beta)/re)* (polystress.tauxz.*((mg.dwdx))+...
%                        polystress.tauyz.*((mg.dwdy)) +...
%                        polystress.tauzz.*(S33F));
%
%
%epFinst=((1-beta)/re)*( tau11F.*conj(S11F)+tau22F.*conj(S22F)+tau33F.*conj(S33F)+...
%		 2*tau12F.*conj(S12F)+2*tau13F.*conj(S13F)+2*tau23F.*conj(S23F));
%
%epinst=((1-beta)/re)*(	(polystress.tauxx).*mg.dudx + (polystress.tauyy).*mg.dvdy + (polystress.tauzz).*mg.dwdz +...
%(polystress.tauxy).*(mg.dudy + mg.dvdx) + (polystress.tauxz).*(mg.dudz +mg.dwdx)+ (polystress.tauyz).*(mg.dwdy+mg.dvdz));
%S11=(mg.dudx);
%S22=(mg.dvdy);
%S33=(mg.dwdz);
%S12=(0.5*( mg.dudy+mg.dvdx));
%S13=(0.5*( mg.dudz+mg.dwdx));
%S23=(0.5*( mg.dvdz+mg.dwdy));
%
%dissinst=(2/re)*(S11.^2+ S22.^2+ S33.^2 + 2*S12.^2+ 2*S13.^2 + 2*S23.^2);

%ep=ep+epinst;
%epF=epF+epFinst;
pixF=pixF+pixFinst;
piyF=piyF+piyFinst;
pizF=pizF+pizFinst;
pix=pix+pixinst;
piy=piy+piyinst;
piz=piz+pizinst;
%dis=dis+dissinst;
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

%dis=dis./nf;
%dis=squeeze(mean(mean(dis,1),2));

fme=fullfile(runFolder,'energy_transferF.mat')
me=matfile(fme,'Writable',true)
%me.ep=ep./nf;
%me.epF=epF./nf;
me.pixF=pixF./nf;
me.piyF=piyF./nf;
me.pizF=pizF./nf;
me.pix=pix./nf;
me.piy=piy./nf;
me.piz=piz./nf;
me.ep=(pix+piy+piz)./nf;
me.epF=(pixF+piyF+pizF)./nf;
%me.dis=dis;
%fny="ygrid.mat";
% % %fmy="/home/skumar67/data-geyink1/skumar67/FV_visco/ygrid.mat";
%fmy=fullfile(runFolder,fny)
%my=matfile(fmy,'Writable',true)
%my.yCheb=yCheb;



