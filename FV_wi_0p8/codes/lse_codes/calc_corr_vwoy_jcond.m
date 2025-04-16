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
phiwu=zeros(Nz,Nx,Ny/2);
phioyu=zeros(Nz,Nx,Ny/2);

phivv=zeros(Nz,Nx,Ny/2);
phiwv=zeros(Nz,Nx,Ny/2);
phioyv=zeros(Nz,Nx,Ny/2);

phivw=zeros(Nz,Nx,Ny/2);
phiww=zeros(Nz,Nx,Ny/2);
phioyw=zeros(Nz,Nx,Ny/2);

phivfx=zeros(Nz,Nx,Ny/2);
phiwfx=zeros(Nz,Nx,Ny/2);
phioyfx=zeros(Nz,Nx,Ny/2);

phivfy=zeros(Nz,Nx,Ny/2);
phiwfy=zeros(Nz,Nx,Ny/2);
phioyfy=zeros(Nz,Nx,Ny/2);

phivfz=zeros(Nz,Nx,Ny/2);
phiwfz=zeros(Nz,Nx,Ny/2);
phioyfz=zeros(Nz,Nx,Ny/2);


vv=0;
vw=0;
voy=0;
ww=0;
woy=0;
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
	
	polyxF=fft2(mt.polyx(:,:,Ny/2+1:end))./(Nz*Nx);
	polyyF=fft2(mt.polyy(:,:,Ny/2+1:end))./(Nz*Nx);
	polyzF=fft2(mt.polyz(:,:,Ny/2+1:end))./(Nz*Nx);

	vfj=m.vFourier(:,:,jcond);
        vfj(1,1)=0;
	wfj=m.wFourier(:,:,jcond);
	wfj(1,1)=0;
	oyfj=fft2(mg.dudz(:,:,jcond)-mg.dwdx(:,:,jcond))./(Nz*Nx);
        oyfj(1,1)=0;

	vv=vv+mean(m.vfield(:,:,jcond).^2,'all');
	vw=vw+mean(m.vfield(:,:,jcond).*m.wfield(:,:,jcond),'all');
        ww=ww+mean(m.wfield(:,:,jcond).*m.wfield(:,:,jcond),'all');

	oy=mg.dudz(:,:,jcond)-mg.dwdx(:,:,jcond)-mean(mg.dudz(:,:,jcond)-mg.dwdx(:,:,jcond),'all' );
	
	voy=voy	+mean(oy.*(m.vfield(:,:,jcond)), 	'all');
	woy=woy	+mean(oy.*(m.wfield(:,:,jcond)), 	'all');
	oyoy=oyoy+mean(oy.^2,'all');
	

	phivu=phivu+conj(vfj).*m.uFourier(:,:,Ny/2+1:end);
	phiwu=phiwu+conj(wfj).*m.uFourier(:,:,Ny/2+1:end);
	phioyu=phioyu+conj(oyfj).*m.uFourier(:,:,Ny/2+1:end);
	
	phivv=phivv+conj(vfj).*m.vFourier(:,:,Ny/2+1:end);
	phiwv=phiwv+conj(wfj).*m.vFourier(:,:,Ny/2+1:end);
        phioyv=phioyv+conj(oyfj).*m.vFourier(:,:,Ny/2+1:end);
	
	phivw=phivw+conj(vfj).*m.wFourier(:,:,Ny/2+1:end);
	phiww=phiww+conj(wfj).*m.wFourier(:,:,Ny/2+1:end);
        phioyw=phioyw+conj(oyfj).*m.wFourier(:,:,Ny/2+1:end);

	phivfx=phivfx+conj(vfj).*polyxF;
	phiwfx=phiwfx+conj(wfj).*polyxF;
        phioyfx=phioyfx+conj(oyfj).*polyxF;

	phivfy=phivfy+conj(vfj).*polyyF;
        phiwfy=phiwfy+conj(wfj).*polyyF;
        phioyfy=phioyfy+conj(oyfj).*polyyF;

	phivfz=phivfz+conj(vfj).*polyzF;
        phiwfz=phiwfz+conj(wfj).*polyzF;
        phioyfz=phioyfz+conj(oyfj).*polyzF;

end

	 phivu= phivu./nf;
         phiwu= phiwu./nf;
        phioyu=phioyu./nf; %+conj(ozfj).*m.uFourier(:,:,Ny/2+1:end);

	 phivv= phivv./nf;
	 phiwv= phiwv./nf;%+conj(vfj).*m.vFourier(:,:,Ny/2+1:end);
        phioyv=phioyv./nf;%+conj(ozfj).*m.vFourier(:,:,Ny/2+1:end);

	 phivw= phivw./nf;
         phiww= phiww./nf;%+conj(vfj).*m.wFourier(:,:,Ny/2+1:end);
        phioyw=phioyw./nf;%	+conj(ozfj).*m.wFourier(:,:,Ny/2+1:end);

	 phivfx= phivfx./nf;
         phiwfx= phiwfx./nf;%+conj(vfj).*polyxF;
        phioyfx=phioyfx./nf;

	 phivfy= phivfy./nf;
         phiwfy= phiwfy./nf;%+conj(vfj).*polyxF;
        phioyfy=phioyfy./nf;

	 phivfz= phivfz./nf;
         phiwfz= phiwfz./nf;%+conj(vfj).*polyxF;
        phioyfz=phioyfz./nf;

	vv=vv./nf;
	vw=vw./nf;
	voy=voy./nf;
	ww=ww./nf;
	woy=woy./nf;
	oyoy=oyoy./nf;
	

Rvu=ifft2(phivu*(Nz*Nx),'symmetric');
Rwu=ifft2(phiwu*(Nz*Nx),'symmetric');
Royu=ifft2(phioyu*(Nz*Nx),'symmetric');

Rvv=ifft2(phivv*(Nz*Nx),'symmetric');
Rwv=ifft2(phiwv*(Nz*Nx),'symmetric');
Royv=ifft2(phioyv*(Nz*Nx),'symmetric');

Rvw=ifft2(phivw*(Nz*Nx),'symmetric');
Rww=ifft2(phiww*(Nz*Nx),'symmetric');
Royw=ifft2(phioyw*(Nz*Nx),'symmetric');

Rvfx=ifft2(phivfx*(Nz*Nx),'symmetric');
Rwfx=ifft2(phiwfx*(Nz*Nx),'symmetric');
Royfx=ifft2(phioyfx*(Nz*Nx),'symmetric');

Rvfy=ifft2(phivfy*(Nz*Nx),'symmetric');
Rwfy=ifft2(phiwfy*(Nz*Nx),'symmetric');
Royfy=ifft2(phioyfy*(Nz*Nx),'symmetric');

Rvfz=ifft2(phivfz*(Nz*Nx),'symmetric');
Rwfz=ifft2(phiwfz*(Nz*Nx),'symmetric');
Royfz=ifft2(phioyfz*(Nz*Nx),'symmetric');

fn=sprintf('vwoy_corr_j_%03d.mat',jcond);
mf=matfile(fn,"Writable",true);

mf.Rvu=Rvu;
mf.Rwu=Rwu;
mf.Royu=Royu;

mf.Rvv=Rvv;
mf.Rwv=Rwv;
mf.Royv=Royv;

mf.Rvw=Rvw;
mf.Rww=Rww;
mf.Royw=Royw;

mf.Rvfx=Rvfx;
mf.Rwfx=Rwfx;
mf.Royfx=Royfx;

mf.Rvfy=Rvfy;
mf.Rwfy=Rwfy;
mf.Royfy=Royfy;

mf.Rvfz=Rvfz;
mf.Rwfz=Rwfz;
mf.Royfz=Royfz;

mf.vv=vv;
mf.vw=vw;
mf.ww=ww;
mf.voy=voy;
mf.woy=woy;
mf.oyoy=oyoy;

mf.j=jc;





