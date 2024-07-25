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
fylse=zeros(nz,nx,ny/2);
fzlse=zeros(nz,nx,ny/2);


uC1=zeros(nz,nx,ny);
vC1=zeros(nz,nx,ny);
wC1=zeros(nz,nx,ny);
fxC1=zeros(nz,nx,ny);
fyC1=zeros(nz,nx,ny);
fzC1=zeros(nz,nx,ny);


uC2=zeros(nz,nx,ny);
vC2=zeros(nz,nx,ny);
wC2=zeros(nz,nx,ny);
fxC2=zeros(nz,nx,ny);
fyC2=zeros(nz,nx,ny);
fzC2=zeros(nz,nx,ny);

lx=4*pi;
lz=2*pi;
xp=lx*[0:nx-1]/nx-lx/2;
zp=lz*[0:nz-1]/nz-lz/2;

load('mean_profiles.mat')
fn=sprintf('vwoy_lse_j_%03d.mat',jcond);
ml=matfile(fn);
%%
v=-0.0263;
w=-0.058;
oy=-1.5;

 ulse=v*ml.L11+w*ml.L12+oy*ml.L13;
 vlse=v*ml.L21+w*ml.L22+oy*ml.L23;
 wlse=v*ml.L31+w*ml.L32+oy*ml.L33;
fxlse=v*ml.L41+w*ml.L42+oy*ml.L43;
fylse=v*ml.L51+w*ml.L52+oy*ml.L53;
fzlse=v*ml.L61+w*ml.L62+oy*ml.L63;

 ulse=(fftshift(fftshift( ulse,1),2));
 vlse=(fftshift(fftshift( vlse,1),2));
 wlse=(fftshift(fftshift( wlse,1),2));
fxlse=(fftshift(fftshift(fxlse,1),2));
fylse=(fftshift(fftshift(fylse,1),2));
fzlse=(fftshift(fftshift(fzlse,1),2));

 uC1(:,:,ny/2+1:end)= ulse;
 vC1(:,:,ny/2+1:end)= vlse;
 wC1(:,:,ny/2+1:end)= wlse;
fxC1(:,:,ny/2+1:end)=fxlse;
fyC1(:,:,ny/2+1:end)=fylse;
fzC1(:,:,ny/2+1:end)=fzlse;

vc1=v;
wc1=w;
oyc1=oy;

%%
v=0.0263;
w=-0.065;
oy=2.3;

 ulse=v*ml.L11+w*ml.L12+oy*ml.L13;
 vlse=v*ml.L21+w*ml.L22+oy*ml.L23;
 wlse=v*ml.L31+w*ml.L32+oy*ml.L33;
fxlse=v*ml.L41+w*ml.L42+oy*ml.L43;
fylse=v*ml.L51+w*ml.L52+oy*ml.L53;
fzlse=v*ml.L61+w*ml.L62+oy*ml.L63;

 ulse=(fftshift(fftshift( ulse,1),2));
 vlse=(fftshift(fftshift( vlse,1),2));
 wlse=(fftshift(fftshift( wlse,1),2));
fxlse=(fftshift(fftshift(fxlse,1),2));
fylse=(fftshift(fftshift(fylse,1),2));
fzlse=(fftshift(fftshift(fzlse,1),2));

 uC2(:,:,ny/2+1:end)= ulse;
 vC2(:,:,ny/2+1:end)= vlse;
 wC2(:,:,ny/2+1:end)= wlse;
fxC2(:,:,ny/2+1:end)=fxlse;
fyC2(:,:,ny/2+1:end)=fylse;
fzC2(:,:,ny/2+1:end)=fzlse;

vc2=v;
wc2=w;
oyc2=oy;
%%
yp=yCheb+1;
[X,Z,Y]=meshgrid(xp,zp,yp);
fn=sprintf('velfield_lse_vwoy_j_%03d.mat',jcond);
m=matfile(fn,'Writable',true);

 m.uQ2= uC1;
 m.vQ2= vC1;
 m.wQ2= wC1;
m.fxQ2=fxC1;
m.fyQ2=fyC1;
m.fzQ2=fzC1;

 m.uQ4= uC2;
 m.vQ4= vC2;
 m.wQ4= wC2;
m.fxQ4=fxC2;
m.fyQ4=fyC2;
m.fzQ4=fzC2;

m.vc1=vc1;
m.wc1=wc1;
m.oyc1=oyc1;

m.vc2=vc2;
m.wc2=wc2;
m.oyc2=oyc2;

m.X=X;
m.Y=Y;
m.Z=Z;
