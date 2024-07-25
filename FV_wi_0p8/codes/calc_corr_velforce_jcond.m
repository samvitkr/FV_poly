close all
clear
load('ygrid.mat')

Nx=512;
Nz=384;
Ny=220;
jcond=156;
jc=jcond-Ny/2;

nf=1;
phiufx=zeros(Nz,Nx,Ny/2);
phivfx=zeros(Nz,Nx,Ny/2);
phiwfx=zeros(Nz,Nx,Ny/2);

tstart=10000;
tend=108000;
tstep=1000;
nf=(tend-tstart)/tstep+1;
%load('lambda_stats.mat')
for time=tstart:tstep:tend
        time
        fvel=sprintf("velfields_%07d.mat",time);
        m=matfile(fvel);
	ft=sprintf("transferfields_%07d.mat",time);
        mt=matfile(ft);
	polyxF=fft2(mt.poly(:,:,Ny/2+1:end))./(Nz*Nx);
%m=matfile('velfields_0050000.mat');
	ufj=m.uFourier(:,:,jcond);
	vfj=m.vFourier(:,:,jcond);
	wfj=m.wFourier(:,:,jcond);
	ufj(1,1)=0;
	vfj(1,1)=0;
	wfj(1,1)=0;
	phiufx=phiufx+conj(ufj).*polyxF;
	phivfx=phivfx+conj(vfj).*polyxF;
	phiwfx=phiwfx+conj(wfj).*polyxF;
end

phiufx=phiufx./nf;
phivfx=phivfx./nf;
phiwfx=phiwfx./nf;

Rufx=ifft2(phiufx*(Nz*Nx),'symmetric');
Rvfx=ifft2(phivfx*(Nz*Nx),'symmetric');
Rwfx=ifft2(phiwfx*(Nz*Nx),'symmetric');

fn=sprintf('vel_force_corr_j_%03d.mat',jcond);
mf=matfile(fn,"Writable",true);
mf.Rufx=Rufx;
mf.Rvfx=Rvfx;
mf.Rwfx=Rwfx;
mf.j=jc;





