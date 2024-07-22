clear;  clc;  close all;
addpath '/home/skumar67/data-geyink1/skumar67/dmsuite'
sourceFolder='../';
% automatically add code modules to matlab path
functionCallStack = dbstack;
[scriptFolder, ~, ~] = fileparts(which(functionCallStack(1).file));
[sourceFolder, ~, ~] = fileparts(scriptFolder);
addpath(fullfile(sourceFolder, 'data_readers'), fullfile(sourceFolder, 'utilities'));
%load('ygrid.mat')
% user input: specify run parameters. This assumes your data is stored in a directory with full path dataFolder/runName
%runName = 'test_wi_6-67_1';
%fileNr = 79000;
nCheb = 220;  % number of Chebyshev collocation points in y (has to be determined a priori)
runFolder='/home/skumar67/data-geyink1/skumar67/FV_wi_1p83'

fy=fullfile(runFolder,'ygrid.mat');
load(fy);
fileNr=70000;
re=4667;
ret=232;
ut=ret/re;
dnu=1/ret;
Nx=512;
Nz=384;
Ny=220;
Lx=  4*pi;
Lz = 2*pi;
xp=Lx*[0:Nx-1]/Nx;
zp=Lz*[0:Nz-1]/Nz;
deltaX=xp(2)-xp(1);
deltaZ=zp(2)-zp(1);
%kx = [0:Nx-1]*pi/Lx;
kx = 2*(pi/Lx)*[0:Nx/2-1, 0, -Nx/2+1:-1];
kz = 2*(pi/Lz)*[0:Nz/2-1, 0, -Nz/2+1:-1];
dkx=kx(2)-kx(1);
dkz=kz(2)-kz(1);
yp= (1-abs(yCheb))./dnu;

[Kx,Kz]=meshgrid(kx,kz);
ft = SpatialFourierTransform();

fng=sprintf("velgrad_%07d.mat",fileNr)
fmg=fullfile(runFolder,fng)
mg=matfile(fmg)

fngt=sprintf("transferfields_%07d.mat",fileNr)
fmgt=fullfile(runFolder,fngt)
mgt=matfile(fmgt)

%fngf=sprintf("velgrad_transfer_flp_%07d.mat",fileNr)
%fmgf=fullfile(runFolder,fngf)
%mgf=matfile(fmgf,'Writable',true)
mgf=matfile('dudxf_test.mat','Writable',true)
mdf=matfile('dragonfly.mat')
D=mdf.D;

[dudxFourier, kx, kz] = ft.transform_to_fourier(mg.dudx, deltaX, deltaZ);
%D=zeros(Nz,Nx,Ny);
%for j=20:95
%	m(j)=exp( -0.2172*log( yp(j) )+0.4452);
%	ka(j)=exp(-1.554*log( yp(j) )+11.44);
%	kb(j)=exp(-0.9912*log(yp(j)) +7.39 );
%	D(:,:,j)=exp(-(  (  abs(Kx)  +m(j)* abs(Kz) ).^2/ka(j)^2  +  (  abs(Kz) -m(j)*abs(Kx) ).^2/kb(j)^2 	));
%	D(:,:,Ny-j+1)=D(:,:,j);
%end


%fdf=fullfile(runFolder,'dragonfly.mat')
%mdf=matfile(fdf);
%mdf.D=D;
%mdf.Kx=Kx;
%mdf.Kz=Kz;
%mdf=matfile('dragonfly.mat');
%D=mdf.D;
mgf.dudx= ft.transform_to_physical( ( mg.dudxF).*D);
mgf.dudx_nof=mg.dudx
%mgf.dvdx= ft.transform_to_physical( ( mg.dvdxF).*D);
%mgf.dwdx= ft.transform_to_physical( ( mg.dwdxF).*D);
%
%mgf.dudy= ft.transform_to_physical( ( mg.dudyF).*D);
%mgf.dvdy= ft.transform_to_physical( ( mg.dvdyF).*D);
%mgf.dwdy= ft.transform_to_physical( ( mg.dwdyF).*D);
%
%mgf.dudz= ft.transform_to_physical( ( mg.dudzF).*D);
%mgf.dvdz= ft.transform_to_physical( ( mg.dvdzF).*D);
%mgf.dwdz= ft.transform_to_physical( ( mg.dwdzF).*D);
%
%mgf.voz=ft.transform_to_physical( ( mgt.vozF).*(D.^2));
%mgf.woy=ft.transform_to_physical( ( mgt.woyF).*(D.^2));
%
%polyF= ft.transform_to_fourier(mgt.poly);
%mgf.poly=ft.transform_to_physical( ( polyF).*(D.^2));
