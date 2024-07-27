close all
clear

Nx=512;
Nz=384;
Ny=220;
jcond=156;
jc=jcond-Ny/2;
ddfilter=zeros(Nz,Nx,Ny/2);
uufilter=zeros(Nz,Nx,Ny/2);
mfil=matfile('../data/filter.mat')
fn=sprintf('../data/vel_corr_j_%03d.mat',jcond);
mf=matfile(fn);
%dfj=mfil.dfil(:,:,jcond);
%ufj=mfil.ufil(:,:,jcond);
%dfj(1,1)=0;
%ufj(1,1)=0;
ddfilter=mfil.dfil(:,:,Ny/2+1:end);
uufilter=mfil.ufil(:,:,Ny/2+1:end);

phiuu=fft2(mf.Ruu);
phivv=fft2(mf.Rvv);
phiww=fft2(mf.Rww);
phiuv=fft2(mf.Ruv);
phivu=fft2(mf.Rvu);
phiuw=fft2(mf.Ruw);
phiwu=fft2(mf.Rwu);
phivw=fft2(mf.Rvw);
phiwv=fft2(mf.Rwv);

fnd=sprintf('../data/vel_corr_ddfilter_j_%03d.mat',jcond);
mfd=matfile(fnd,'Writable',true);

mfd.Ruu=ifft2(ddfilter.*phiuu,'symmetric');
mfd.Rvv=ifft2(ddfilter.*phivv,'symmetric');
mfd.Rww=ifft2(ddfilter.*phiww,'symmetric');
mfd.Ruv=ifft2(ddfilter.*phiuv,'symmetric');
mfd.Rvu=ifft2(ddfilter.*phivu,'symmetric');
mfd.Ruw=ifft2(ddfilter.*phiuw,'symmetric');
mfd.Rwu=ifft2(ddfilter.*phiwu,'symmetric');
mfd.Rvw=ifft2(ddfilter.*phivw,'symmetric');
mfd.Rwv=ifft2(ddfilter.*phiwv,'symmetric');


fnu=sprintf('../data/vel_corr_uufilter_j_%03d.mat',jcond);
mfu=matfile(fnu,'Writable',true);

mfu.Ruu=ifft2(uufilter.*phiuu,'symmetric');
mfu.Rvv=ifft2(uufilter.*phivv,'symmetric');
mfu.Rww=ifft2(uufilter.*phiww,'symmetric');
mfu.Ruv=ifft2(uufilter.*phiuv,'symmetric');
mfu.Rvu=ifft2(uufilter.*phivu,'symmetric');
mfu.Ruw=ifft2(uufilter.*phiuw,'symmetric');
mfu.Rwu=ifft2(uufilter.*phiwu,'symmetric');
mfu.Rvw=ifft2(uufilter.*phivw,'symmetric');
mfu.Rwv=ifft2(uufilter.*phiwv,'symmetric');
