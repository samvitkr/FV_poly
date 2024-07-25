close all
clear
load('ygrid.mat')

Nx=512;
Nz=384;
Ny=220;
jcond=188;
jc=jcond-Ny/2;

nf=1;

phivu=zeros(Nz,Nx,Ny/2);
phiozu=zeros(Nz,Nx,Ny/2);

phivv=zeros(Nz,Nx,Ny/2);
phiozv=zeros(Nz,Nx,Ny/2);

phivw=zeros(Nz,Nx,Ny/2);
phiozw=zeros(Nz,Nx,Ny/2);

phivfx=zeros(Nz,Nx,Ny/2);
phiozfx=zeros(Nz,Nx,Ny/2);

ozoz=0;

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
	vfj=m.vFourier(:,:,jcond);
	vfj(1,1)=0;
	ozfj=fft2(mg.dvdx(:,:,jcond)-mg.dudy(:,:,jcond))./(Nz*Nx);
        ozfj(1,1)=0;

	oz=mg.dvdx(:,:,jcond)-mg.dudy(:,:,jcond) - mean( mg.dvdx(:,:,jcond)-mg.dudy(:,:,jcond),'all' );
	ozoz=ozoz+mean(oz.^2,'all');

	phivfx=phivfx+conj(vfj).*polyxF;
	phiozfx=phiozfx+conj(ozfj).*polyxF;
	
	phivu=phivu+conj(vfj).*m.uFourier(:,:,Ny/2+1:end);
	phiozu=phiozu+conj(ozfj).*m.uFourier(:,:,Ny/2+1:end);
	
	phivv=phivv+conj(vfj).*m.vFourier(:,:,Ny/2+1:end);
        phiozv=phiozv+conj(ozfj).*m.vFourier(:,:,Ny/2+1:end);
	
	phivw=phivw+conj(vfj).*m.wFourier(:,:,Ny/2+1:end);
        phiozw=phiozw+conj(ozfj).*m.wFourier(:,:,Ny/2+1:end);
end

	phivfx=phivfx./nf;%+conj(vfj).*polyxF;
        phiozfx=phiozfx./nf;

        phivu=phivu./nf;
        phiozu=phiozu./nf; %+conj(ozfj).*m.uFourier(:,:,Ny/2+1:end);

        phivv=phivv./nf;%+conj(vfj).*m.vFourier(:,:,Ny/2+1:end);
        phiozv=phiozv./nf;%+conj(ozfj).*m.vFourier(:,:,Ny/2+1:end);

        phivw=phivw./nf;%+conj(vfj).*m.wFourier(:,:,Ny/2+1:end);
        phiozw=phiozw./nf;%	+conj(ozfj).*m.wFourier(:,:,Ny/2+1:end);

	ozoz=ozoz./nf;

Rvfx=ifft2(phivfx*(Nz*Nx),'symmetric');
Rozfx=ifft2(phiozfx*(Nz*Nx),'symmetric');

Rvu=ifft2(phivu*(Nz*Nx),'symmetric');
Rozu=ifft2(phiozu*(Nz*Nx),'symmetric');

Rvv=ifft2(phivv*(Nz*Nx),'symmetric');
Rozv=ifft2(phiozv*(Nz*Nx),'symmetric');

Rvw=ifft2(phivw*(Nz*Nx),'symmetric');
Rozw=ifft2(phiozw*(Nz*Nx),'symmetric');

fn=sprintf('voz_woy_corr_j_%03d.mat',jcond);
mf=matfile(fn,"Writable",true);
mf.Rvfx=Rvfx;
mf.Rozfx=Rozfx;

mf.Rvu=Rvu;
mf.Rozu=Rozu;

mf.Rvv=Rvv;
mf.Rozv=Rozv;

mf.Rvw=Rvw;
mf.Rozw=Rozw;

mf.ozoz=ozoz;
mf.j=jc;





