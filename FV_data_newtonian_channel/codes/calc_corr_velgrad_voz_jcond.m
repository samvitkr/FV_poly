close all
clear
%load('ygrid.mat')

Nx=512;
Nz=384;
Ny=220;
jcond=156;
jc=jcond-Ny/2;

nf=1;

phivu=zeros(Nz,Nx,Ny/2);
phiozu=zeros(Nz,Nx,Ny/2);

phivv=zeros(Nz,Nx,Ny/2);
phiozv=zeros(Nz,Nx,Ny/2);

phivw=zeros(Nz,Nx,Ny/2);
phiozw=zeros(Nz,Nx,Ny/2);

phivdudx=zeros(Nz,Nx,Ny/2);
phivdvdx=zeros(Nz,Nx,Ny/2);
phivdwdx=zeros(Nz,Nx,Ny/2);

phivdudy=zeros(Nz,Nx,Ny/2);
phivdvdy=zeros(Nz,Nx,Ny/2);
phivdwdy=zeros(Nz,Nx,Ny/2);

phivdudz=zeros(Nz,Nx,Ny/2);
phivdvdz=zeros(Nz,Nx,Ny/2);
phivdwdz=zeros(Nz,Nx,Ny/2);

phiozdudx=zeros(Nz,Nx,Ny/2);
phiozdvdx=zeros(Nz,Nx,Ny/2);
phiozdwdx=zeros(Nz,Nx,Ny/2);

phiozdudy=zeros(Nz,Nx,Ny/2);
phiozdvdy=zeros(Nz,Nx,Ny/2);
phiozdwdy=zeros(Nz,Nx,Ny/2);

phiozdudz=zeros(Nz,Nx,Ny/2);
phiozdvdz=zeros(Nz,Nx,Ny/2);
phiozdwdz=zeros(Nz,Nx,Ny/2);

ozoz=0;

tstart=0000;
tend=100000;
tstep=1000;
nf=(tend-tstart)/tstep+1;
%load('lambda_stats.mat')
for time=tstart:tstep:tend
        time
        fvel=sprintf("../data/velfields_%07d.mat",time);
        m=matfile(fvel);
	ft=sprintf("../data/transferfields_%07d.mat",time);
        mt=matfile(ft);
	fg=sprintf("../data/velgrad_%07d.mat",time);
        mg=matfile(fg);
	
%	polyxF=fft2(mt.poly(:,:,Ny/2+1:end))./(Nz*Nx);
	vfj=m.vFourier(:,:,jcond);
	vfj(1,1)=0;
	ozfj=fft2(mg.dvdx(:,:,jcond)-mg.dudy(:,:,jcond))./(Nz*Nx);
        ozfj(1,1)=0;

	oz=mg.dvdx(:,:,jcond)-mg.dudy(:,:,jcond) - mean( mg.dvdx(:,:,jcond)-mg.dudy(:,:,jcond),'all' );
	ozoz=ozoz+mean(oz.^2,'all');

%	phivfx=phivfx+conj(vfj).*polyxF;
%	phiozfx=phiozfx+conj(ozfj).*polyxF;
	phivu=phivu+conj(vfj).*m.uFourier(:,:,Ny/2+1:end);
	phiozu=phiozu+conj(ozfj).*m.uFourier(:,:,Ny/2+1:end);
	phivv=phivv+conj(vfj).*m.vFourier(:,:,Ny/2+1:end);
        phiozv=phiozv+conj(ozfj).*m.vFourier(:,:,Ny/2+1:end);
	phivw=phivw+conj(vfj).*m.wFourier(:,:,Ny/2+1:end);
        phiozw=phiozw+conj(ozfj).*m.wFourier(:,:,Ny/2+1:end);

	phivdudx=phivdudx+conj(vfj).*fft2(mg.dudx(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdvdx=phivdvdx+conj(vfj).*fft2(mg.dvdx(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdwdx=phivdwdx+conj(vfj).*fft2(mg.dwdx(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdudy=phivdudy+conj(vfj).*fft2(mg.dudy(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdvdy=phivdvdy+conj(vfj).*fft2(mg.dvdy(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdwdy=phivdwdy+conj(vfj).*fft2(mg.dwdy(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdudz=phivdudz+conj(vfj).*fft2(mg.dudz(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdvdz=phivdvdz+conj(vfj).*fft2(mg.dvdz(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdwdz=phivdwdz+conj(vfj).*fft2(mg.dwdz(:,:,Ny/2+1:end))./(Nx*Nz);

%	phiozdudx=phiozdudx+conj(ozfj).*fft2(mg.dudx(:,:,Ny/2+1:end))./(Nx*Nz);
        phiozdvdx=phiozdvdx+conj(ozfj).*fft2(mg.dvdx(:,:,Ny/2+1:end))./(Nx*Nz);
        phiozdwdx=phiozdwdx+conj(ozfj).*fft2(mg.dwdx(:,:,Ny/2+1:end))./(Nx*Nz);
        phiozdudy=phiozdudy+conj(ozfj).*fft2(mg.dudy(:,:,Ny/2+1:end))./(Nx*Nz);
        phiozdvdy=phiozdvdy+conj(ozfj).*fft2(mg.dvdy(:,:,Ny/2+1:end))./(Nx*Nz);
        phiozdwdy=phiozdwdy+conj(ozfj).*fft2(mg.dwdy(:,:,Ny/2+1:end))./(Nx*Nz);
        phiozdudz=phiozdudz+conj(ozfj).*fft2(mg.dudz(:,:,Ny/2+1:end))./(Nx*Nz);
        phiozdvdz=phiozdvdz+conj(ozfj).*fft2(mg.dvdz(:,:,Ny/2+1:end))./(Nx*Nz);
        phiozdwdz=phiozdwdz+conj(ozfj).*fft2(mg.dwdz(:,:,Ny/2+1:end))./(Nx*Nz);

end

%	phivfx=phivfx./nf;%+conj(vfj).*polyxF;
%        phiozfx=phiozfx./nf;

        phivu=phivu./nf;
        phiozu=phiozu./nf; %+conj(ozfj).*m.uFourier(:,:,Ny/2+1:end);

        phivv=phivv./nf;%+conj(vfj).*m.vFourier(:,:,Ny/2+1:end);
        phiozv=phiozv./nf;%+conj(ozfj).*m.vFourier(:,:,Ny/2+1:end);

        phivw=phivw./nf;%+conj(vfj).*m.wFourier(:,:,Ny/2+1:end);
        phiozw=phiozw./nf;%	+conj(ozfj).*m.wFourier(:,:,Ny/2+1:end);

	ozoz=ozoz./nf;

	phivdudx=phivdudx./nf;
	phivdvdx=phivdvdx./nf;
	phivdwdx=phivdwdx./nf;
	phivdudy=phivdudy./nf;
	phivdvdy=phivdvdy./nf;
	phivdwdy=phivdwdy./nf;
	phivdudz=phivdudz./nf;
	phivdvdz=phivdvdz./nf;
	phivdwdz=phivdwdz./nf;
	
	phiozdudx=phiozdudx./nf;
	phiozdvdx=phiozdvdx./nf;
	phiozdwdx=phiozdwdx./nf;
	phiozdudy=phiozdudy./nf;
	phiozdvdy=phiozdvdy./nf;
	phiozdwdy=phiozdwdy./nf;
	phiozdudz=phiozdudz./nf;
	phiozdvdz=phiozdvdz./nf;
	phiozdwdz=phiozdwdz./nf;


%Rvfx=ifft2(phivfx*(Nz*Nx),'symmetric');
%Rozfx=ifft2(phiozfx*(Nz*Nx),'symmetric');

Rvu=ifft2(phivu*(Nz*Nx),'symmetric');
Rozu=ifft2(phiozu*(Nz*Nx),'symmetric');

Rvv=ifft2(phivv*(Nz*Nx),'symmetric');
Rozv=ifft2(phiozv*(Nz*Nx),'symmetric');

Rvw=ifft2(phivw*(Nz*Nx),'symmetric');
Rozw=ifft2(phiozw*(Nz*Nx),'symmetric');

Rvdudx=ifft2(phivdudx*(Nz*Nx),'symmetric');
Rvdvdx=ifft2(phivdvdx*(Nz*Nx),'symmetric');
Rvdwdx=ifft2(phivdwdx*(Nz*Nx),'symmetric');
Rvdudy=ifft2(phivdudy*(Nz*Nx),'symmetric');
Rvdvdy=ifft2(phivdvdy*(Nz*Nx),'symmetric');
Rvdwdy=ifft2(phivdwdy*(Nz*Nx),'symmetric');
Rvdudz=ifft2(phivdudz*(Nz*Nx),'symmetric');
Rvdvdz=ifft2(phivdvdz*(Nz*Nx),'symmetric');
Rvdwdz=ifft2(phivdwdz*(Nz*Nx),'symmetric');

Rozdudx=ifft2(phiozdudx*(Nz*Nx),'symmetric');
Rozdvdx=ifft2(phiozdvdx*(Nz*Nx),'symmetric');
Rozdwdx=ifft2(phiozdwdx*(Nz*Nx),'symmetric');
Rozdudy=ifft2(phiozdudy*(Nz*Nx),'symmetric');
Rozdvdy=ifft2(phiozdvdy*(Nz*Nx),'symmetric');
Rozdwdy=ifft2(phiozdwdy*(Nz*Nx),'symmetric');
Rozdudz=ifft2(phiozdudz*(Nz*Nx),'symmetric');
Rozdvdz=ifft2(phiozdvdz*(Nz*Nx),'symmetric');
Rozdwdz=ifft2(phiozdwdz*(Nz*Nx),'symmetric');

fn=sprintf('../data/voz_velgrad_corr_j_%03d.mat',jcond);
mf=matfile(fn,"Writable",true);
%mf.Rvfx=Rvfx;
%mf.Rozfx=Rozfx;

mf.Rvu=Rvu;
mf.Rozu=Rozu;

mf.Rvv=Rvv;
mf.Rozv=Rozv;

mf.Rvw=Rvw;
mf.Rozw=Rozw;

mf.ozoz=ozoz;
mf.j=jc;

mf.Rozdudx=Rozdudx;
mf.Rozdvdx=Rozdvdx;
mf.Rozdwdx=Rozdwdx;
mf.Rozdudy=Rozdudy;
mf.Rozdvdy=Rozdvdy;
mf.Rozdwdy=Rozdwdy;
mf.Rozdudz=Rozdudz;
mf.Rozdvdz=Rozdvdz;
mf.Rozdwdz=Rozdwdz;

mf.Rvdudx=Rvdudx;
mf.Rvdvdx=Rvdvdx;
mf.Rvdwdx=Rvdwdx;
mf.Rvdudy=Rvdudy;
mf.Rvdvdy=Rvdvdy;
mf.Rvdwdy=Rvdwdy;
mf.Rvdudz=Rvdudz;
mf.Rvdvdz=Rvdvdz;
mf.Rvdwdz=Rvdwdz;
