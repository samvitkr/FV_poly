close all
clear
load('../data/ygrid.mat')

Nx=512;
Nz=384;
Ny=220;
jcond=156;
jc=Ny-jcond+1;

nf=1;
phivv=zeros(Nz,Nx,Ny/2);
phivu=zeros(Nz,Nx,Ny/2);
phivw=zeros(Nz,Nx,Ny/2);

phivdudx=zeros(Nz,Nx,Ny/2);
phivdvdx=zeros(Nz,Nx,Ny/2);
phivdwdx=zeros(Nz,Nx,Ny/2);

phivdudy=zeros(Nz,Nx,Ny/2);
phivdvdy=zeros(Nz,Nx,Ny/2);
phivdwdy=zeros(Nz,Nx,Ny/2);

phivdudz=zeros(Nz,Nx,Ny/2);
phivdvdz=zeros(Nz,Nx,Ny/2);
phivdwdz=zeros(Nz,Nx,Ny/2);

phivfx=zeros(Nz,Nx,Ny/2);
phivfy=zeros(Nz,Nx,Ny/2);
phivfz=zeros(Nz,Nx,Ny/2);
phivvoz=zeros(Nz,Nx,Ny/2);
phivwoy=zeros(Nz,Nx,Ny/2);

tstart=10000;
tend=108000;
tstep=1000;
nf=2*((tend-tstart-tstep)/tstep+1);
%load('lambda_stats.mat')
for time=tstart+tstep:tstep:tend
        time

        fvel=sprintf("../data/velfields_%07d.mat",time-tstep);
        m=matfile(fvel);
	vfj=m.vFourier(:,:,jcond);
        vfj(1,1)=0;
        vfjt=m.vFourier(:,:,jc);
        vfjt(1,1)=0;
	clear m

	fvel=sprintf("../data/velfields_%07d.mat",time);
        m=matfile(fvel);
        phivu=phivu+conj(vfj).*m.uFourier(:,:,Ny/2+1:end)-flip(conj(vfjt).*m.uFourier(:,:,1:Ny/2),3);
        phivv=phivv+conj(vfj).*m.vFourier(:,:,Ny/2+1:end)+flip(conj(vfjt).*m.vFourier(:,:,1:Ny/2),3);        
	phivw=phivw+conj(vfj).*m.wFourier(:,:,Ny/2+1:end)-flip(conj(vfjt).*m.wFourier(:,:,1:Ny/2),3);
	clear m

	ft=sprintf("../data/transferfields_%07d.mat",time);
        mt=matfile(ft);

        polyxF=       fft2(mt.polyx(:,:,Ny/2+1:end))./(Nz*Nx);
	polyxFt=flip( fft2(mt.polyx(:,:,    1:Ny/2))./(Nz*Nx) ,3);
	polyyF=       fft2(mt.polyy(:,:,Ny/2+1:end))./(Nz*Nx);
        polyyFt=flip( fft2(mt.polyy(:,:,    1:Ny/2))./(Nz*Nx) ,3);
	polyzF=       fft2(mt.polyz(:,:,Ny/2+1:end))./(Nz*Nx);
        polyzFt=flip( fft2(mt.polyz(:,:,    1:Ny/2))./(Nz*Nx) ,3);

	vozF=       fft2(mt.voz(:,:,Ny/2+1:end))./(Nz*Nx);
        vozFt=flip( fft2(mt.voz(:,:,    1:Ny/2))./(Nz*Nx) ,3);
	woyF=       fft2(mt.woy(:,:,Ny/2+1:end))./(Nz*Nx);
        woyFt=flip( fft2(mt.woy(:,:,    1:Ny/2))./(Nz*Nx) ,3);

	phivfx=phivfx+conj(vfj).*polyxF-conj(vfjt).*polyxFt;
	phivfy=phivfy+conj(vfj).*polyyF+conj(vfjt).*polyyFt;
	phivfz=phivfz+conj(vfj).*polyzF-conj(vfjt).*polyzFt;

	phivvoz=phivvoz+conj(vfj).*vozF-conj(vfjt).*vozFt;
	phivwoy=phivwoy+conj(vfj).*woyF-conj(vfjt).*woyFt;
	clear mt

	fvelg=sprintf("../data/velgrad_%07d.mat",time);
        mg=matfile(fvelg);

        phivdudx=phivdudx+conj(vfj).*fft2(mg.dudx(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdvdx=phivdvdx+conj(vfj).*fft2(mg.dvdx(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdwdx=phivdwdx+conj(vfj).*fft2(mg.dwdx(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdudy=phivdudy+conj(vfj).*fft2(mg.dudy(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdvdy=phivdvdy+conj(vfj).*fft2(mg.dvdy(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdwdy=phivdwdy+conj(vfj).*fft2(mg.dwdy(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdudz=phivdudz+conj(vfj).*fft2(mg.dudz(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdvdz=phivdvdz+conj(vfj).*fft2(mg.dvdz(:,:,Ny/2+1:end))./(Nx*Nz);
        phivdwdz=phivdwdz+conj(vfj).*fft2(mg.dwdz(:,:,Ny/2+1:end))./(Nx*Nz);

	phivdudx=phivdudx-conj(vfjt).*flip(fft2(mg.dudx(:,:,1:Ny/2))./(Nx*Nz),3);
        phivdvdx=phivdvdx+conj(vfjt).*flip(fft2(mg.dvdx(:,:,1:Ny/2))./(Nx*Nz),3);
        phivdwdx=phivdwdx-conj(vfjt).*flip(fft2(mg.dwdx(:,:,1:Ny/2))./(Nx*Nz),3);

	phivdudy=phivdudy+conj(vfjt).*flip(fft2(mg.dudy(:,:,1:Ny/2))./(Nx*Nz),3);
        phivdvdy=phivdvdy-conj(vfjt).*flip(fft2(mg.dvdy(:,:,1:Ny/2))./(Nx*Nz),3);
        phivdwdy=phivdwdy+conj(vfjt).*flip(fft2(mg.dwdy(:,:,1:Ny/2))./(Nx*Nz),3);

	phivdudz=phivdudz-conj(vfjt).*flip(fft2(mg.dudz(:,:,1:Ny/2))./(Nx*Nz),3);
        phivdvdz=phivdvdz+conj(vfjt).*flip(fft2(mg.dvdz(:,:,1:Ny/2))./(Nx*Nz),3);
        phivdwdz=phivdwdz-conj(vfjt).*flip(fft2(mg.dwdz(:,:,1:Ny/2))./(Nx*Nz),3);

	clear mg
end

phivv=phivv./nf;
phivu=phivu./nf;
phivw=phivw./nf;

phivfx=phivfx./nf;
phivfy=phivfy./nf;
phivfz=phivfz./nf;
phivvoz=phivvoz./nf;
phivwoy=phivwoy./nf;

phivdudx=phivdudx./nf;
phivdvdx=phivdvdx./nf;
phivdwdx=phivdwdx./nf;

phivdudy=phivdudy./nf;
phivdvdy=phivdvdy./nf;
phivdwdy=phivdwdy./nf;

phivdudz=phivdudz./nf;
phivdvdz=phivdvdz./nf;
phivdwdz=phivdwdz./nf;

fn=sprintf('../data/vel_corr_reflect_F_j_%03d.mat',jcond);
mf=matfile(fn,"Writable",true);

mf.Rvv=ifft2(phivv*(Nz*Nx),'symmetric');clear phivv
mf.Rvu=ifft2(phivu*(Nz*Nx),'symmetric');clear phivu
mf.Rvw=ifft2(phivw*(Nz*Nx),'symmetric');clear phivw

mf.Rvfx=ifft2(phivfx*(Nz*Nx),'symmetric');clear phivfx
mf.Rvfy=ifft2(phivfy*(Nz*Nx),'symmetric');clear phivfy
mf.Rvfz=ifft2(phivfz*(Nz*Nx),'symmetric');clear phivfz

mf.Rvvoz=ifft2(phivvoz*(Nz*Nx),'symmetric');clear phivvoz
mf.Rvwoy=ifft2(phivwoy*(Nz*Nx),'symmetric');clear phivwoy

mf.Rvdudx=ifft2(phivdudx*(Nz*Nx),'symmetric');clear phivdudx
mf.Rvdvdx=ifft2(phivdvdx*(Nz*Nx),'symmetric');clear phivdvdx
mf.Rvdwdx=ifft2(phivdwdx*(Nz*Nx),'symmetric');clear phivdwdx

mf.Rvdudy=ifft2(phivdudy*(Nz*Nx),'symmetric');clear phivdudy
mf.Rvdvdy=ifft2(phivdvdy*(Nz*Nx),'symmetric');clear phivdvdy
mf.Rvdwdy=ifft2(phivdwdy*(Nz*Nx),'symmetric');clear phivdwdy

mf.Rvdudz=ifft2(phivdudz*(Nz*Nx),'symmetric');clear phivdudz
mf.Rvdvdz=ifft2(phivdvdz*(Nz*Nx),'symmetric');clear phivdvdz
mf.Rvdwdz=ifft2(phivdwdz*(Nz*Nx),'symmetric');clear phivdwdz
