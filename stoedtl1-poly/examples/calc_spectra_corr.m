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
Wi=0.8;
beta=0.9;
L_max=100;
a=1-3/L_max^2;
re=4667
runFolder='/home/skumar67/data-geyink1/skumar67/FV_wi_0p8'
%runFolder='/home/skumar67/data-geyink1/skumar67/FV_data_'
[xGridPointsDns, yGridPointsDns, ~, ~] = read_grid(runFolder, ngx, ngy);
deltaX = xGridPointsDns(2) - xGridPointsDns(1);
deltaZ = Lz / (ngz - 1);
[yCheb, ~] = chebdif(nCheb, 1);
interpY = InterpolationOperatorY(yGridPointsDns, yCheb);
ft = SpatialFourierTransform();
ep=zeros(ngz-1,ngx-1,nCheb);
epF=zeros(ngz-1,ngx-1,nCheb);
phiuu=ep;
phivv=ep;
phiww=ep;
phiuv=ep;
phiuw=ep;
phivw=ep;

%Ruu=ep;
%Rvv=ep;
%Rww=ep;
%Ruv=ep;
%Ruw=ep;
%Rvw=ep;

fstart=10000;
fend=108000;
fskip=1000;
nf=(fend-fstart)/fskip+1;
for fileNr=fstart:fskip:fend
fng=sprintf("velfields_%07d.mat",fileNr)
fmg=fullfile(runFolder,fng)
mg=matfile(fmg)
phiuu=phiuu+mg.uFourier.*conj(mg.uFourier);
phivv=phivv+mg.vFourier.*conj(mg.vFourier);
phiww=phiww+mg.wFourier.*conj(mg.wFourier);

phiuv=phiuv+mg.uFourier.*conj(mg.vFourier);
phiuw=phiuw+mg.uFourier.*conj(mg.wFourier);
phivw=phivw+mg.vFourier.*conj(mg.wFourier);
end
fme=fullfile(runFolder,'spectra_corr.mat')
me=matfile(fme,'Writable',true)
me.phiuu=phiuu./nf;
me.phivv=phivv./nf;
me.phiww=phiww./nf;
me.phiuv=phiuv./nf;
me.phiuw=phiuw./nf;
me.phivw=phivw./nf;

me.Ruu=ft.transform_to_physical( phiuu );
me.Rvv=ft.transform_to_physical( phivv );
me.Rww=ft.transform_to_physical( phiww );
me.Ruv=ft.transform_to_physical( phiuv );
me.Ruw=ft.transform_to_physical( phiuw );
me.Rvw=ft.transform_to_physical( phivw );

