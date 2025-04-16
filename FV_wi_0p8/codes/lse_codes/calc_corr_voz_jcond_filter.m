close all
clear
%load('ygrid.mat')
%
Nx=512;
Nz=384;
Ny=220;
jcond=188;
%jc=jcond-Ny/2;
%
%nf=1;
%
%phivu=zeros(Nz,Nx,Ny/2);
%phiozu=zeros(Nz,Nx,Ny/2);
phiozoz=zeros(Nz,Nx);
%phivv=zeros(Nz,Nx,Ny/2);
%phiozv=zeros(Nz,Nx,Ny/2);
%
%phivw=zeros(Nz,Nx,Ny/2);
%phiozw=zeros(Nz,Nx,Ny/2);
%
%phivfx=zeros(Nz,Nx,Ny/2);
%phiozfx=zeros(Nz,Nx,Ny/2);
%
ozoz=0;
%
tstart=10000;
tend=108000;
tstep=1000;
nf=(tend-tstart)/tstep+1;
% %load('lambda_stats.mat')
for time=tstart:tstep:tend
        time
        fvel=sprintf("../data/velfields_%07d.mat",time);
        m=matfile(fvel);
%	ft=sprintf("transferfields_%07d.mat",time);
%        mt=matfile(ft);
	fg=sprintf("../data/velgrad_%07d.mat",time);
        mg=matfile(fg);
%	
%	polyxF=fft2(mt.poly(:,:,Ny/2+1:end))./(Nz*Nx);
%	vfj=m.vFourier(:,:,jcond);
%	vfj(1,1)=0;
	ozfj=fft2(mg.dvdx(:,:,jcond)-mg.dudy(:,:,jcond))./(Nz*Nx);
	phiozoz=phiozoz+conj(ozfj).*(ozfj);
%        ozfj(1,1)=0;
%
%	oz=mg.dvdx(:,:,jcond)-mg.dudy(:,:,jcond) - mean( mg.dvdx(:,:,jcond)-mg.dudy(:,:,jcond),'all' );
%	ozoz=ozoz+mean(oz.^2,'all');
%
%	phivfx=phivfx+conj(vfj).*polyxF;
%	phiozfx=phiozfx+conj(ozfj).*polyxF;
%	
%	phivu=phivu+conj(vfj).*m.uFourier(:,:,Ny/2+1:end);
%	phiozu=phiozu+conj(ozfj).*m.uFourier(:,:,Ny/2+1:end);
%	
%	phivv=phivv+conj(vfj).*m.vFourier(:,:,Ny/2+1:end);
%        phiozv=phiozv+conj(ozfj).*m.vFourier(:,:,Ny/2+1:end);
%	
%	phivw=phivw+conj(vfj).*m.wFourier(:,:,Ny/2+1:end);
%        phiozw=phiozw+conj(ozfj).*m.wFourier(:,:,Ny/2+1:end);
end
phiozoz=phiozoz./nf;
%
%	phivfx=phivfx./nf;%+conj(vfj).*polyxF;
%        phiozfx=phiozfx./nf;
%
%        phivu=phivu./nf;
%        phiozu=phiozu./nf; %+conj(ozfj).*m.uFourier(:,:,Ny/2+1:end);
%
%        phivv=phivv./nf;%+conj(vfj).*m.vFourier(:,:,Ny/2+1:end);
%        phiozv=phiozv./nf;%+conj(ozfj).*m.vFourier(:,:,Ny/2+1:end);
%
%        phivw=phivw./nf;%+conj(vfj).*m.wFourier(:,:,Ny/2+1:end);
%        phiozw=phiozw./nf;%	+conj(ozfj).*m.wFourier(:,:,Ny/2+1:end);
%
%	ozoz=ozoz./nf;
%
%Rvfx=ifft2(phivfx*(Nz*Nx),'symmetric');
%Rozfx=ifft2(phiozfx*(Nz*Nx),'symmetric');
%
%Rvu=ifft2(phivu*(Nz*Nx),'symmetric');
%Rozu=ifft2(phiozu*(Nz*Nx),'symmetric');
%
%Rvv=ifft2(phivv*(Nz*Nx),'symmetric');
%Rozv=ifft2(phiozv*(Nz*Nx),'symmetric');
%
%Rvw=ifft2(phivw*(Nz*Nx),'symmetric');
%Rozw=ifft2(phiozw*(Nz*Nx),'symmetric');
%
%fn=sprintf('voz_woy_corr_j_%03d.mat',jcond);
%mf=matfile(fn,"Writable",true);
%mf.Rvfx=Rvfx;
%mf.Rozfx=Rozfx;
%
%mf.Rvu=Rvu;
%mf.Rozu=Rozu;
%
%mf.Rvv=Rvv;
%mf.Rozv=Rozv;
%
%mf.Rvw=Rvw;
%mf.Rozw=Rozw;
%
%mf.ozoz=ozoz;
%mf.j=jc;
%close all
%clear
load('../data/ygrid.mat')
%Nx=512;
%Nz=384;
%Ny=220;
%jcond=188;
jc=jcond-Ny/2;
ddfilter=zeros(Nz,Nx,Ny/2);
uufilter=zeros(Nz,Nx,Ny/2);
mfil=matfile('../data/filter.mat')
%fn=sprintf('../data/voz_woy_corr_j_%03d.mat',jcond);
fn=sprintf('../data/voz_corr_j_%03d.mat',jcond);
mf=matfile(fn)
%dfj=mfil.dfil(:,:,jcond);
%ufj=mfil.ufil(:,:,jcond);
%dfj(1,1)=0;
%ufj(1,1)=0;
ddfilter=mfil.dfil(:,:,Ny/2+1:end);
uufilter=mfil.ufil(:,:,Ny/2+1:end);
ozozdd=sum(phiozoz.*ddfilter(jc),'all');
ozozuu=sum(phiozoz.*uufilter(jc),'all');

ozoz=sum(phiozoz,'all');
vv=mf.Rvv(1,1,jc);
ozv=mf.Rozv(1,1,jc);

phivv=fft2(mf.Rvv);
phivu=fft2(mf.Rvu);
phivw=fft2(mf.Rvw);
phivfx=fft2(mf.Rvfx);
phiozv=fft2(mf.Rozv);
phiozu=fft2(mf.Rozu);
phiozw=fft2(mf.Rozw);
phiozfx=fft2(mf.Rozfx);
fnd=sprintf('../data/voz_corr_ddfilter_j_%03d.mat',jcond);
mfd=matfile(fnd,'Writable',true);

mfd.Rvv=ifft2(ddfilter.*phivv,'symmetric');
mfd.Rvu=ifft2(ddfilter.*phivu,'symmetric');
mfd.Rvw=ifft2(ddfilter.*phivw,'symmetric');
mfd.Rvfx=ifft2(ddfilter.*phivfx,'symmetric');

mfd.Rozv=ifft2(ddfilter.*phiozv,'symmetric');
mfd.Rozu=ifft2(ddfilter.*phiozu,'symmetric');
mfd.Rozw=ifft2(ddfilter.*phiozw,'symmetric');
mfd.Rozfx=ifft2(ddfilter.*phiozfx,'symmetric');
mfd.ozozdd=ozozdd;

mfd.ozoz=ozoz;
mfd.vv=vv;
mfd.ozv=ozv;

mfd.yCheb=yCheb(Ny/2+1:end);
mfd.j=jc;

fnu=sprintf('../data/voz_corr_uufilter_j_%03d.mat',jcond);
mfu=matfile(fnu,'Writable',true);

mfu.Rvv=ifft2(uufilter.*phivv,'symmetric');
mfu.Rvu=ifft2(uufilter.*phivu,'symmetric');
mfu.Rvw=ifft2(uufilter.*phivw,'symmetric');
mfu.Rvfx=ifft2(uufilter.*phivfx,'symmetric');

mfu.Rozv=ifft2(uufilter.*phiozv,'symmetric');
mfu.Rozu=ifft2(uufilter.*phiozu,'symmetric');
mfu.Rozw=ifft2(uufilter.*phiozw,'symmetric');
mfu.Rozfx=ifft2(uufilter.*phiozfx,'symmetric');

mfu.ozozuu=ozozuu;

mfu.ozoz=ozoz;
mfu.vv=vv;
mfu.ozv=ozv;

mfu.yCheb=yCheb(Ny/2+1:end);
mfu.j=jc;
