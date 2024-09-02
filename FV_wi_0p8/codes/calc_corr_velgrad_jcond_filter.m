close all
clear
load('../data/ygrid.mat')
Nx=512;
Nz=384;
Ny=220;
jcond=156;
jc=jcond-Ny/2;
ddfilter=zeros(Nz,Nx,Ny/2);
uufilter=zeros(Nz,Nx,Ny/2);
mfil=matfile('../data/filter.mat')
fn=sprintf('../data/velgrad_corr_j_%03d.mat',jcond);
mf=matfile(fn);
%dfj=mfil.dfil(:,:,jcond);
%ufj=mfil.ufil(:,:,jcond);
%dfj(1,1)=0;
%ufj(1,1)=0;
ddfilter=mfil.dfil(:,:,Ny/2+1:end);
uufilter=mfil.ufil(:,:,Ny/2+1:end);

phiuu=fft2(mf.Ruu);
phivv=fft2(mf.Rvv);
phiuv=fft2(mf.Ruv);
phivu=fft2(mf.Rvu);
phiuw=fft2(mf.Ruw);
phivw=fft2(mf.Rvw);

phiududx=fft2(mf.Rududx);
phiudvdx=fft2(mf.Rudvdx);
phiudwdx=fft2(mf.Rudwdx);
phiududy=fft2(mf.Rududy);
phiudvdy=fft2(mf.Rudvdy);
phiudwdy=fft2(mf.Rudwdy);
phiududz=fft2(mf.Rududz);
phiudvdz=fft2(mf.Rudvdz);
phiudwdz=fft2(mf.Rudwdz);

phivdudx=fft2(mf.Rvdudx);
phivdvdx=fft2(mf.Rvdvdx);
phivdwdx=fft2(mf.Rvdwdx);
phivdudy=fft2(mf.Rvdudy);
phivdvdy=fft2(mf.Rvdvdy);
phivdwdy=fft2(mf.Rvdwdy);
phivdudz=fft2(mf.Rvdudz);
phivdvdz=fft2(mf.Rvdvdz);
phivdwdz=fft2(mf.Rvdwdz);


fnd=sprintf('../data/velgrad_corr_ddfilter_j_%03d.mat',jcond);
mfd=matfile(fnd,'Writable',true);

mfd.Ruu=ifft2(ddfilter.*phiuu,'symmetric');
mfd.Rvv=ifft2(ddfilter.*phivv,'symmetric');
mfd.Ruv=ifft2(ddfilter.*phiuv,'symmetric');
mfd.Rvu=ifft2(ddfilter.*phivu,'symmetric');
mfd.Ruw=ifft2(ddfilter.*phiuw,'symmetric');
mfd.Rvw=ifft2(ddfilter.*phivw,'symmetric');

mfd.Rududx=ifft2(ddfilter.*phiududx,'symmetric');
mfd.Rudvdx=ifft2(ddfilter.*phiudvdx,'symmetric');
mfd.Rudwdx=ifft2(ddfilter.*phiudwdx,'symmetric');
mfd.Rududy=ifft2(ddfilter.*phiududy,'symmetric');
mfd.Rudvdy=ifft2(ddfilter.*phiudvdy,'symmetric');
mfd.Rudwdy=ifft2(ddfilter.*phiudwdy,'symmetric');
mfd.Rududz=ifft2(ddfilter.*phiududz,'symmetric');
mfd.Rudvdz=ifft2(ddfilter.*phiudvdz,'symmetric');
mfd.Rudwdz=ifft2(ddfilter.*phiudwdz,'symmetric');

mfd.Rvdudx=ifft2(ddfilter.*phivdudx,'symmetric');
mfd.Rvdvdx=ifft2(ddfilter.*phivdvdx,'symmetric');
mfd.Rvdwdx=ifft2(ddfilter.*phivdwdx,'symmetric');
mfd.Rvdudy=ifft2(ddfilter.*phivdudy,'symmetric');
mfd.Rvdvdy=ifft2(ddfilter.*phivdvdy,'symmetric');
mfd.Rvdwdy=ifft2(ddfilter.*phivdwdy,'symmetric');
mfd.Rvdudz=ifft2(ddfilter.*phivdudz,'symmetric');
mfd.Rvdvdz=ifft2(ddfilter.*phivdvdz,'symmetric');
mfd.Rvdwdz=ifft2(ddfilter.*phivdwdz,'symmetric');

mfd.yCheb=yCheb(Ny/2+1:end);
mfd.j=jc;

fnu=sprintf('../data/velgrad_corr_uufilter_j_%03d.mat',jcond);
mfu=matfile(fnu,'Writable',true);

mfu.Ruu=ifft2(uufilter.*phiuu,'symmetric');
mfu.Rvv=ifft2(uufilter.*phivv,'symmetric');
mfu.Ruv=ifft2(uufilter.*phiuv,'symmetric');
mfu.Rvu=ifft2(uufilter.*phivu,'symmetric');
mfu.Ruw=ifft2(uufilter.*phiuw,'symmetric');
mfu.Rvw=ifft2(uufilter.*phivw,'symmetric');

mfu.Rududx=ifft2(uufilter.*phiududx,'symmetric');
mfu.Rudvdx=ifft2(uufilter.*phiudvdx,'symmetric');
mfu.Rudwdx=ifft2(uufilter.*phiudwdx,'symmetric');
mfu.Rududy=ifft2(uufilter.*phiududy,'symmetric');
mfu.Rudvdy=ifft2(uufilter.*phiudvdy,'symmetric');
mfu.Rudwdy=ifft2(uufilter.*phiudwdy,'symmetric');
mfu.Rududz=ifft2(uufilter.*phiududz,'symmetric');
mfu.Rudvdz=ifft2(uufilter.*phiudvdz,'symmetric');
mfu.Rudwdz=ifft2(uufilter.*phiudwdz,'symmetric');

mfu.Rvdudx=ifft2(uufilter.*phivdudx,'symmetric');
mfu.Rvdvdx=ifft2(uufilter.*phivdvdx,'symmetric');
mfu.Rvdwdx=ifft2(uufilter.*phivdwdx,'symmetric');
mfu.Rvdudy=ifft2(uufilter.*phivdudy,'symmetric');
mfu.Rvdvdy=ifft2(uufilter.*phivdvdy,'symmetric');
mfu.Rvdwdy=ifft2(uufilter.*phivdwdy,'symmetric');
mfu.Rvdudz=ifft2(uufilter.*phivdudz,'symmetric');
mfu.Rvdvdz=ifft2(uufilter.*phivdvdz,'symmetric');
mfu.Rvdwdz=ifft2(uufilter.*phivdwdz,'symmetric');

mfu.yCheb=yCheb(Ny/2+1:end);
mfu.j=jc;
