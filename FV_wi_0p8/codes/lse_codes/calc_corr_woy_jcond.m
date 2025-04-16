close all
clear
load('ygrid.mat')

Nx=512;
Nz=384;
Ny=220;
jcond=188;
jc=jcond-Ny/2;

nf=1;

phiwu=zeros(Nz,Nx,Ny/2);
phioyu=zeros(Nz,Nx,Ny/2);

phiwv=zeros(Nz,Nx,Ny/2);
phioyv=zeros(Nz,Nx,Ny/2);

phiww=zeros(Nz,Nx,Ny/2);
phioyw=zeros(Nz,Nx,Ny/2);

phiwfx=zeros(Nz,Nx,Ny/2);
phioyfx=zeros(Nz,Nx,Ny/2);

oyoy=0;

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
	fg=sprintf("velgrad_%07d.mat",time);
        mg=matfile(fg);
	
	polyxF=fft2(mt.poly(:,:,Ny/2+1:end))./(Nz*Nx);
	wfj=m.wFourier(:,:,jcond);
	wfj(1,1)=0;
	oyfj=fft2(mg.dudz(:,:,jcond)-mg.dwdx(:,:,jcond))./(Nz*Nx);
        oyfj(1,1)=0;

	oy=mg.dudz(:,:,jcond)-mg.dwdx(:,:,jcond) - mean( mg.dudz(:,:,jcond)-mg.dwdx(:,:,jcond),'all' );
	oyoy=oyoy+mean(oy.^2,'all');

	phiwfx=phiwfx+conj(wfj).*polyxF;
	phioyfx=phioyfx+conj(oyfj).*polyxF;
	
	phiwu=phiwu+conj(wfj).*m.uFourier(:,:,Ny/2+1:end);
	phioyu=phioyu+conj(oyfj).*m.uFourier(:,:,Ny/2+1:end);
	
	phiwv=phiwv+conj(wfj).*m.vFourier(:,:,Ny/2+1:end);
        phioyv=phioyv+conj(oyfj).*m.vFourier(:,:,Ny/2+1:end);
	
	phiww=phiww+conj(wfj).*m.wFourier(:,:,Ny/2+1:end);
        phioyw=phioyw+conj(oyfj).*m.wFourier(:,:,Ny/2+1:end);
end

	phiwfx=phiwfx./nf;%+conj(vfj).*polyxF;
        phioyfx=phioyfx./nf;

        phiwu=phiwu./nf;
        phioyu=phioyu./nf; %+conj(ozfj).*m.uFourier(:,:,Ny/2+1:end);

        phiwv=phiwv./nf;%+conj(vfj).*m.vFourier(:,:,Ny/2+1:end);
        phioyv=phioyv./nf;%+conj(ozfj).*m.vFourier(:,:,Ny/2+1:end);

        phiww=phiww./nf;%+conj(vfj).*m.wFourier(:,:,Ny/2+1:end);
        phioyw=phioyw./nf;%	+conj(ozfj).*m.wFourier(:,:,Ny/2+1:end);

	oyoy=oyoy./nf;

Rwfx=ifft2(phiwfx*(Nz*Nx),'symmetric');
Royfx=ifft2(phioyfx*(Nz*Nx),'symmetric');

Rwu=ifft2(phiwu*(Nz*Nx),'symmetric');
Royu=ifft2(phioyu*(Nz*Nx),'symmetric');

Rwv=ifft2(phiwv*(Nz*Nx),'symmetric');
Royv=ifft2(phioyv*(Nz*Nx),'symmetric');

Rww=ifft2(phiww*(Nz*Nx),'symmetric');
Royw=ifft2(phioyw*(Nz*Nx),'symmetric');

fn=sprintf('voz_woy_corr_j_%03d.mat',jcond);
mf=matfile(fn,"Writable",true);
mf.Rwfx=Rwfx;
mf.Royfx=Royfx;

mf.Rwu=Rwu;
mf.Royu=Royu;

mf.Rwv=Rwv;
mf.Royv=Royv;

mf.Rww=Rww;
mf.Royw=Royw;

mf.oyoy=oyoy;
mf.j=jc;





