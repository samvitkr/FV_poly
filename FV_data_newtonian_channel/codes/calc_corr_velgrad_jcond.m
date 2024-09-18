close all
clear
load('../data/ygrid.mat')

Nx=512;
Nz=384;
Ny=220;
jcond=156;
jc=jcond-Ny/2;

nf=1;
phiuu=zeros(Nz,Nx,Ny/2);
phivv=zeros(Nz,Nx,Ny/2);
%phiww=zeros(Nz,Nx,Ny/2);
phiuv=zeros(Nz,Nx,Ny/2);
phivu=zeros(Nz,Nx,Ny/2);
phiuw=zeros(Nz,Nx,Ny/2);
%phiwu=zeros(Nz,Nx,Ny/2);
phivw=zeros(Nz,Nx,Ny/2);
%phiwv=zeros(Nz,Nx,Ny/2);

phiududx=zeros(Nz,Nx,Ny/2);
phiudvdx=zeros(Nz,Nx,Ny/2);
phiudwdx=zeros(Nz,Nx,Ny/2);

phiududy=zeros(Nz,Nx,Ny/2);
phiudvdy=zeros(Nz,Nx,Ny/2);
phiudwdy=zeros(Nz,Nx,Ny/2);

phiududz=zeros(Nz,Nx,Ny/2);
phiudvdz=zeros(Nz,Nx,Ny/2);
phiudwdz=zeros(Nz,Nx,Ny/2);

phivdudx=zeros(Nz,Nx,Ny/2);
phivdvdx=zeros(Nz,Nx,Ny/2);
phivdwdx=zeros(Nz,Nx,Ny/2);

phivdudy=zeros(Nz,Nx,Ny/2);
phivdvdy=zeros(Nz,Nx,Ny/2);
phivdwdy=zeros(Nz,Nx,Ny/2);

phivdudz=zeros(Nz,Nx,Ny/2);
phivdvdz=zeros(Nz,Nx,Ny/2);
phivdwdz=zeros(Nz,Nx,Ny/2);


tstart=0000;
tend=100000;
tstep=1000;
nf=(tend-tstart)/tstep+1;
%load('lambda_stats.mat')
for time=tstart:tstep:tend
        time
        fvel=sprintf("../data/velfields_%07d.mat",time);
        m=matfile(fvel);
	fvelg=sprintf("../data/velgrad_%07d.mat",time);
        mg=matfile(fvelg);

%m=matfile('velfields_0050000.mat');
	ufj=m.uFourier(:,:,jcond);
	vfj=m.vFourier(:,:,jcond);
%	wfj=m.wFourier(:,:,jcond);
	ufj(1,1)=0;
	vfj(1,1)=0;
%	wfj(1,1)=0;

	phiuu=phiuu+conj(ufj).*m.uFourier(:,:,Ny/2+1:end);
	phivv=phivv+conj(vfj).*m.vFourier(:,:,Ny/2+1:end);
%	phiww=phiww+conj(wfj).*m.wFourier(:,:,Ny/2+1:end);
	phiuv=phiuv+conj(ufj).*m.vFourier(:,:,Ny/2+1:end);
	phivu=phivu+conj(vfj).*m.uFourier(:,:,Ny/2+1:end);
	phiuw=phiuw+conj(ufj).*m.wFourier(:,:,Ny/2+1:end);
%	phiwu=phiwu+conj(wfj).*m.uFourier(:,:,Ny/2+1:end);
	phivw=phivw+conj(vfj).*m.wFourier(:,:,Ny/2+1:end);
%	phiwv=phiwv+conj(wfj).*m.vFourier(:,:,Ny/2+1:end);

        phiududx=phiududx+conj(ufj).*fft2(mg.dudx(:,:,Ny/2+1:end))./(Nx*Nz);
	phiudvdx=phiudvdx+conj(ufj).*fft2(mg.dvdx(:,:,Ny/2+1:end))./(Nx*Nz);
	phiudwdx=phiudwdx+conj(ufj).*fft2(mg.dwdx(:,:,Ny/2+1:end))./(Nx*Nz);
	phiududy=phiududy+conj(ufj).*fft2(mg.dudy(:,:,Ny/2+1:end))./(Nx*Nz);
	phiudvdy=phiudvdy+conj(ufj).*fft2(mg.dvdy(:,:,Ny/2+1:end))./(Nx*Nz);
	phiudwdy=phiudwdy+conj(ufj).*fft2(mg.dwdy(:,:,Ny/2+1:end))./(Nx*Nz);
	phiududz=phiududz+conj(ufj).*fft2(mg.dudz(:,:,Ny/2+1:end))./(Nx*Nz);
	phiudvdz=phiudvdz+conj(ufj).*fft2(mg.dvdz(:,:,Ny/2+1:end))./(Nx*Nz);
	phiudwdz=phiudwdz+conj(ufj).*fft2(mg.dwdz(:,:,Ny/2+1:end))./(Nx*Nz);	

        phivdudx=phivdudx+conj(vfj).*fft2(mg.dudx(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdvdx=phivdvdx+conj(vfj).*fft2(mg.dvdx(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdwdx=phivdwdx+conj(vfj).*fft2(mg.dwdx(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdudy=phivdudy+conj(vfj).*fft2(mg.dudy(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdvdy=phivdvdy+conj(vfj).*fft2(mg.dvdy(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdwdy=phivdwdy+conj(vfj).*fft2(mg.dwdy(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdudz=phivdudz+conj(vfj).*fft2(mg.dudz(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdvdz=phivdvdz+conj(vfj).*fft2(mg.dvdz(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdwdz=phivdwdz+conj(vfj).*fft2(mg.dwdz(:,:,Ny/2+1:end))./(Nx*Nz);

end

phiuu=phiuu./nf;
phivv=phivv./nf;
%phiww=phiww./nf;
phiuv=phiuv./nf;
phivu=phivu./nf;
phiuw=phiuw./nf;
%phiwu=phiwu./nf;
phivw=phivw./nf;
%phiwv=phiwv./nf;

phiududx=phiududx./nf;
phiudvdx=phiudvdx./nf;
phiudwdx=phiudwdx./nf;
phiududy=phiududy./nf;
phiudvdy=phiudvdy./nf;
phiudwdy=phiudwdy./nf;
phiududz=phiududz./nf;
phiudvdz=phiudvdz./nf;
phiudwdz=phiudwdz./nf;

phivdudx=phivdudx./nf;
phivdvdx=phivdvdx./nf;
phivdwdx=phivdwdx./nf;
phivdudy=phivdudy./nf;
phivdvdy=phivdvdy./nf;
phivdwdy=phivdwdy./nf;
phivdudz=phivdudz./nf;
phivdvdz=phivdvdz./nf;
phivdwdz=phivdwdz./nf;


Ruu=ifft2(phiuu*(Nz*Nx),'symmetric');
Rvv=ifft2(phivv*(Nz*Nx),'symmetric');
%Rww=ifft2(phiww*(Nz*Nx),'symmetric');
Ruv=ifft2(phiuv*(Nz*Nx),'symmetric');
Rvu=ifft2(phivu*(Nz*Nx),'symmetric');
Ruw=ifft2(phiuw*(Nz*Nx),'symmetric');
%Rwu=ifft2(phiwu*(Nz*Nx),'symmetric');
Rvw=ifft2(phivw*(Nz*Nx),'symmetric');
%Rwv=ifft2(phiwv*(Nz*Nx),'symmetric');

Rududx=ifft2(phiududx*(Nz*Nx),'symmetric');
Rudvdx=ifft2(phiudvdx*(Nz*Nx),'symmetric');
Rudwdx=ifft2(phiudwdx*(Nz*Nx),'symmetric');
Rududy=ifft2(phiududy*(Nz*Nx),'symmetric');
Rudvdy=ifft2(phiudvdy*(Nz*Nx),'symmetric');
Rudwdy=ifft2(phiudwdy*(Nz*Nx),'symmetric');
Rududz=ifft2(phiududz*(Nz*Nx),'symmetric');
Rudvdz=ifft2(phiudvdz*(Nz*Nx),'symmetric');
Rudwdz=ifft2(phiudwdz*(Nz*Nx),'symmetric');

Rvdudx=ifft2(phivdudx*(Nz*Nx),'symmetric');
Rvdvdx=ifft2(phivdvdx*(Nz*Nx),'symmetric');
Rvdwdx=ifft2(phivdwdx*(Nz*Nx),'symmetric');
Rvdudy=ifft2(phivdudy*(Nz*Nx),'symmetric');
Rvdvdy=ifft2(phivdvdy*(Nz*Nx),'symmetric');
Rvdwdy=ifft2(phivdwdy*(Nz*Nx),'symmetric');
Rvdudz=ifft2(phivdudz*(Nz*Nx),'symmetric');
Rvdvdz=ifft2(phivdvdz*(Nz*Nx),'symmetric');
Rvdwdz=ifft2(phivdwdz*(Nz*Nx),'symmetric');


fn=sprintf('../data/velgrad_corr_j_%03d.mat',jcond);
mf=matfile(fn,"Writable",true);
mf.Ruu=Ruu;
mf.Rvv=Rvv;
%mf.Rww=Rww;
mf.Ruv=Ruv;
mf.Rvu=Rvu;
mf.Ruw=Ruw;
%mf.Rwu=Rwu;
mf.Rvw=Rvw;
%mf.Rwv=Rwv;
mf.yCheb=yCheb(Ny/2+1:end);
mf.j=jc;

mf.Rududx=Rududx;
mf.Rudvdx=Rudvdx;
mf.Rudwdx=Rudwdx;
mf.Rududy=Rududy;
mf.Rudvdy=Rudvdy;
mf.Rudwdy=Rudwdy;
mf.Rududz=Rududz;
mf.Rudvdz=Rudvdz;
mf.Rudwdz=Rudwdz;

mf.Rvdudx=Rvdudx;
mf.Rvdvdx=Rvdvdx;
mf.Rvdwdx=Rvdwdx;
mf.Rvdudy=Rvdudy;
mf.Rvdvdy=Rvdvdy;
mf.Rvdwdy=Rvdwdy;
mf.Rvdudz=Rvdudz;
mf.Rvdvdz=Rvdvdz;
mf.Rvdwdz=Rvdwdz;


