close all
clear
jcond=188;

nx=512;
nz=384;
ny=220;
ulse=zeros(nz,nx,ny/2);
vlse=zeros(nz,nx,ny/2);
wlse=zeros(nz,nx,ny/2);
fxlse=zeros(nz,nx,ny/2);

uC1=zeros(nz,nx,ny);
vC1=zeros(nz,nx,ny);
wC1=zeros(nz,nx,ny);
fxC1=zeros(nz,nx,ny);

uC2=zeros(nz,nx,ny);
vC2=zeros(nz,nx,ny);
wC2=zeros(nz,nx,ny);
fxC2=zeros(nz,nx,ny);

lx=4*pi;
lz=2*pi;
xp=lx*[0:nx-1]/nx-lx/2;
zp=lz*[0:nz-1]/nz-lz/2;

load('mean_profiles.mat')
fn=sprintf('woy_lse_j_%03d.mat',jcond);
ml=matfile(fn);
%%
w=-0.08;
oy=-1.62;

 ulse=w*ml.L11+oy*ml.L12;%+w*ml.L13;
 vlse=w*ml.L21+oy*ml.L22;%+w*ml.L23;
 wlse=w*ml.L31+oy*ml.L32;%+w*ml.L33;
fxlse=w*ml.L41+oy*ml.L42;%+w*ml.L43;

ulse=(fftshift(fftshift(ulse,1),2));
vlse=(fftshift(fftshift(vlse,1),2));
wlse=(fftshift(fftshift(wlse,1),2));
fxlse=(fftshift(fftshift(fxlse,1),2));
uC1(:,:,ny/2+1:end)=ulse;
vC1(:,:,ny/2+1:end)=vlse;
wC1(:,:,ny/2+1:end)=wlse;
fxC1(:,:,ny/2+1:end)=fxlse;
%%
w=0.055;
oy=-1.66;

 ulse=w*ml.L11+oy*ml.L12;%+w*ml.L13;
 vlse=w*ml.L21+oy*ml.L22;%+w*ml.L23;
 wlse=w*ml.L31+oy*ml.L32;%+w*ml.L33;
fxlse=w*ml.L41+oy*ml.L42;%+w*ml.L43;

ulse=(fftshift(fftshift(ulse,1),2));
vlse=(fftshift(fftshift(vlse,1),2));
wlse=(fftshift(fftshift(wlse,1),2));
fxlse=(fftshift(fftshift(fxlse,1),2));
uC2(:,:,ny/2+1:end)=ulse;
vC2(:,:,ny/2+1:end)=vlse;
wC2(:,:,ny/2+1:end)=wlse;
fxC2(:,:,ny/2+1:end)=fxlse;
%%
yp=yCheb+1;
[X,Z,Y]=meshgrid(xp,zp,yp);
fn=sprintf('velfield_lse_woy_j_%03d.mat',jcond);
m=matfile(fn,'Writable',true);
m.uQ2=uC1;
m.vQ2=vC1;
m.wQ2=wC1;
m.uQ4=uC2;
m.vQ4=vC2;
m.wQ4=wC2;
m.X=X;
m.Y=Y;
m.Z=Z;
m.fxQ2=fxC1;
m.fxQ4=fxC2;
