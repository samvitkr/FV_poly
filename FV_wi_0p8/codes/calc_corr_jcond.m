close all
clear
load('ygrid.mat')

Nx=512;
Nz=384;
Ny=220;
jcond=188;
jc=jcond-Ny/2;

nf=1;
phiuu=zeros(Nz,Nx,Ny/2);
phivv=zeros(Nz,Nx,Ny/2);
phiww=zeros(Nz,Nx,Ny/2);

phiuv=zeros(Nz,Nx,Ny/2);
phivu=zeros(Nz,Nx,Ny/2);

phiuw=zeros(Nz,Nx,Ny/2);
phiwu=zeros(Nz,Nx,Ny/2);

phivw=zeros(Nz,Nx,Ny/2);
phiwv=zeros(Nz,Nx,Ny/2);

tstart=10000;
tend=108000;
tstep=1000;
nf=(tend-tstart)/tstep+1;
%load('lambda_stats.mat')
for time=tstart:tstep:tend
        time
        fvel=sprintf("velfields_%07d.mat",time);
        m=matfile(fvel);

%m=matfile('velfields_0050000.mat');
	ufj=m.uFourier(:,:,jcond);
	vfj=m.vFourier(:,:,jcond);
	wfj=m.wFourier(:,:,jcond);
	ufj(1,1)=0;
	vfj(1,1)=0;
	wfj(1,1)=0;

	phiuu=phiuu+conj(ufj).*m.uFourier(:,:,Ny/2+1:end);
	phivv=phivv+conj(vfj).*m.vFourier(:,:,Ny/2+1:end);
	phiww=phiww+conj(wfj).*m.wFourier(:,:,Ny/2+1:end);
	
	phiuv=phiuv+conj(ufj).*m.vFourier(:,:,Ny/2+1:end);
	phivu=phivu+conj(vfj).*m.uFourier(:,:,Ny/2+1:end);
	
	phiuw=phiuw+conj(ufj).*m.wFourier(:,:,Ny/2+1:end);
	phiwu=phiwu+conj(wfj).*m.uFourier(:,:,Ny/2+1:end);
	
	phivw=phivw+conj(vfj).*m.wFourier(:,:,Ny/2+1:end);
	phiwv=phiwv+conj(wfj).*m.vFourier(:,:,Ny/2+1:end);
end

phiuu=phiuu./nf;
phivv=phivv./nf;
phiww=phiww./nf;

phiuv=phiuv./nf;
phivu=phivu./nf;

phiuw=phiuw./nf;
phiwu=phiwu./nf;

phivw=phivw./nf;
phiwv=phiwv./nf;

Ruu=ifft2(phiuu*(Nz*Nx),'symmetric');
Rvv=ifft2(phivv*(Nz*Nx),'symmetric');
Rww=ifft2(phiww*(Nz*Nx),'symmetric');

Ruv=ifft2(phiuv*(Nz*Nx),'symmetric');
Rvu=ifft2(phivu*(Nz*Nx),'symmetric');

Ruw=ifft2(phiuw*(Nz*Nx),'symmetric');
Rwu=ifft2(phiwu*(Nz*Nx),'symmetric');

Rvw=ifft2(phivw*(Nz*Nx),'symmetric');
Rwv=ifft2(phiwv*(Nz*Nx),'symmetric');

fn=sprintf('vel_corr_j_%03d.mat',jcond);
mf=matfile(fn,"Writable",true);
mf.Ruu=Ruu;
mf.Rvv=Rvv;
mf.Rww=Rww;

mf.Ruv=Ruv;
mf.Rvu=Rvu;

mf.Ruw=Ruw;
mf.Rwu=Rwu;

mf.Rvw=Rvw;
mf.Rwv=Rwv;
mf.yCheb=yCheb(Ny/2+1:end);
mf.j=jc;





