close all
clear
jcond=156;

nx=512;
nz=384;
ny=220;
ulse=zeros(nz,nx,ny/2);
vlse=zeros(nz,nx,ny/2);
wlse=zeros(nz,nx,ny/2);
fxlse=zeros(nz,nx,ny/2);

uQ2=zeros(nz,nx,ny);
vQ2=zeros(nz,nx,ny);
wQ2=zeros(nz,nx,ny);
fxQ2=zeros(nz,nx,ny);

uQ4=zeros(nz,nx,ny);
vQ4=zeros(nz,nx,ny);
wQ4=zeros(nz,nx,ny);
fxQ4=zeros(nz,nx,ny);

lx=4*pi;
lz=2*pi;
xp=lx*[0:nx-1]/nx-lx/2;
zp=lz*[0:nz-1]/nz-lz/2;

load('mean_profiles.mat')
fn=sprintf('vel_lse_j_%03d.mat',jcond);
ml=matfile(fn);
%%
u=-0.115;
v=0.065;
w=0;
ulse=u*ml.L11+v*ml.L12+w*ml.L13;
vlse=u*ml.L21+v*ml.L22+w*ml.L23;
wlse=u*ml.L31+v*ml.L32+w*ml.L33;
fxlse=u*ml.L41+v*ml.L42+w*ml.L43;

ulse=(fftshift(fftshift(ulse,1),2));
vlse=(fftshift(fftshift(vlse,1),2));
wlse=(fftshift(fftshift(wlse,1),2));
fxlse=(fftshift(fftshift(fxlse,1),2));
uQ2(:,:,ny/2+1:end)=ulse;
vQ2(:,:,ny/2+1:end)=vlse;
wQ2(:,:,ny/2+1:end)=wlse;
fxQ2(:,:,ny/2+1:end)=fxlse;
%%
u=0.095;
v=-0.05;
w=0;
ulse=u*ml.L11+v*ml.L12+w*ml.L13;
vlse=u*ml.L21+v*ml.L22+w*ml.L23;
wlse=u*ml.L31+v*ml.L32+w*ml.L33;
fxlse=u*ml.L41+v*ml.L42+w*ml.L43;
ulse=(fftshift(fftshift(ulse,1),2));
vlse=(fftshift(fftshift(vlse,1),2));
wlse=(fftshift(fftshift(wlse,1),2));
fxlse=(fftshift(fftshift(fxlse,1),2));
uQ4(:,:,ny/2+1:end)=ulse;
vQ4(:,:,ny/2+1:end)=vlse;
wQ4(:,:,ny/2+1:end)=wlse;
fxQ4(:,:,ny/2+1:end)=fxlse;
%%
yp=yCheb+1;
[X,Z,Y]=meshgrid(xp,zp,yp);
fn=sprintf('velfield_lse_j_%03d.mat',jcond);
m=matfile(fn,'Writable',true);
m.uQ2=uQ2;
m.vQ2=vQ2;
m.wQ2=wQ2;
m.uQ4=uQ4;
m.vQ4=vQ4;
m.wQ4=wQ4;
m.X=X;
m.Y=Y;
m.Z=Z;
m.fxQ2=fxQ2;
m.fxQ4=fxQ4;


