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
runFolder='/home/skumar67/data-geyink1/skumar67/FV_wi_13p5'

fy=fullfile(runFolder,'ygrid.mat');
load(fy);
fileNr=79000;
re=4667;
ret=154;
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


fnm=sprintf("velfields_%07d.mat",fileNr);
fm=fullfile(runFolder,fnm)
m=matfile(fm,'Writable',true)

%fng=sprintf("velgrad_%07d.mat",fileNr)
%fmg=fullfile(runFolder,fng)
%mg=matfile(fmg)
%
%fngt=sprintf("transferfields_%07d.mat",fileNr)
%fmgt=fullfile(runFolder,fngt)
%mgt=matfile(fmgt)

fngf=sprintf("velfield_flp_%07d.mat",fileNr)
fmgf=fullfile(runFolder,fngf)
mgf=matfile(fmgf,'Writable',true)

fmdname=fullfile(runFolder,'dragonfly.mat')
md=matfile(fmdname)
%md=matfile('dragonfly.mat')
D=md.D;

%[dudxF, kx, kz] = ft.transform_to_fourier(mg.dudx, deltaX, deltaZ);
%[dvdxF] = ft.transform_to_fourier(mg.dvdx);
%[dwdxF] = ft.transform_to_fourier(mg.dwdx);
%
%[dudyF] = ft.transform_to_fourier(mg.dudy);
%[dvdyF] = ft.transform_to_fourier(mg.dvdy);
%[dwdyF] = ft.transform_to_fourier(mg.dwdy);
%
%[dudzF] = ft.transform_to_fourier(mg.dudz);
%[dvdzF] = ft.transform_to_fourier(mg.dvdz);
%[dwdzF] = ft.transform_to_fourier(mg.dwdz);

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
mgf.u= single(ft.transform_to_physical( ( m.uFourier).*D));
mgf.v= single(ft.transform_to_physical( ( m.vFourier).*D));
mgf.w= single(ft.transform_to_physical( ( m.wFourier).*D));

%mgf.dudy= single(ft.transform_to_physical( ( dudyF).*D));
%mgf.dvdy= single(ft.transform_to_physical( ( dvdyF).*D));
%mgf.dwdy= single(ft.transform_to_physical( ( dwdyF).*D));
%
%mgf.dudz= single(ft.transform_to_physical( ( dudzF).*D));
%mgf.dvdz= single(ft.transform_to_physical( ( dvdzF).*D));
%mgf.dwdz= single(ft.transform_to_physical( ( dwdzF).*D));
%
%mgf.voz=single(ft.transform_to_physical( ( mgt.vozF).*(D.^2)));
%mgf.woy=single(ft.transform_to_physical( ( mgt.woyF).*(D.^2)));
%
%polyF= ft.transform_to_fourier(mgt.poly);
%mgf.poly=ft.transform_to_physical( ( polyF).*(D.^2));
