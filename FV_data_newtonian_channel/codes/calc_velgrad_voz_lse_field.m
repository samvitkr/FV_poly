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

load('../data/mean_profiles.mat')
fn=sprintf('../data/velgrad_voz_lse_j_%03d.mat',jcond);
ml=matfile(fn);
%%
%jcond=156
v=0.0482;
oz=-0.6863;

   ulse=fftshift(fftshift(v*ml.L11+  oz*ml.L12,1),2);...+w*ml.L13;
   vlse=fftshift(fftshift(v*ml.L21+  oz*ml.L22,1),2);...+w*ml.L23;
   wlse=fftshift(fftshift(v*ml.L31+  oz*ml.L32,1),2);...+w*ml.L33;

dudxlse=fftshift(fftshift(v*ml.L41+  oz*ml.L42,1),2);
dvdxlse=fftshift(fftshift(v*ml.L51+  oz*ml.L52,1),2);
dwdxlse=fftshift(fftshift(v*ml.L61+  oz*ml.L62,1),2);

dudylse=fftshift(fftshift(v*ml.L71+  oz*ml.L72,1),2);
dvdylse=fftshift(fftshift(v*ml.L81+  oz*ml.L82,1),2);
dwdylse=fftshift(fftshift(v*ml.L91+  oz*ml.L92,1),2);

dudzlse=fftshift(fftshift(v*ml.L101 +oz*ml.L102,1),2);
dvdzlse=fftshift(fftshift(v*ml.L111 +oz*ml.L112,1),2);
dwdzlse=fftshift(fftshift(v*ml.L121 +oz*ml.L122,1),2);


yp=yCheb+1;
[X,Z,Y]=meshgrid(xp,zp,yp(ny/2+1:end));
fn=sprintf('../data/velgrad_voz_field_lseQ2_j_%03d.mat',jcond);
m=matfile(fn,'Writable',true);
m.ozond=oz;
m.vcond=v;
m.u=ulse;
m.v=vlse;
m.w=wlse;
m.dudx=dudxlse;
m.dvdx=dvdxlse;
m.dwdx=dwdxlse;
m.dudy=dudylse;
m.dvdy=dvdylse;
m.dwdy=dwdylse;
m.dudz=dudzlse;
m.dvdz=dvdzlse;
m.dwdz=dwdzlse;
m.X=X;
m.Y=Y;
m.Z=Z;

%%
%jcond=156
v=-0.0514;
oz=-0.4667;

   ulse=fftshift(fftshift(v*ml.L11+  oz*ml.L12,1),2);...+w*ml.L13;
   vlse=fftshift(fftshift(v*ml.L21+  oz*ml.L22,1),2);...+w*ml.L23;
   wlse=fftshift(fftshift(v*ml.L31+  oz*ml.L32,1),2);...+w*ml.L33;

dudxlse=fftshift(fftshift(v*ml.L41+  oz*ml.L42,1),2);
dvdxlse=fftshift(fftshift(v*ml.L51+  oz*ml.L52,1),2);
dwdxlse=fftshift(fftshift(v*ml.L61+  oz*ml.L62,1),2);

dudylse=fftshift(fftshift(v*ml.L71+  oz*ml.L72,1),2);
dvdylse=fftshift(fftshift(v*ml.L81+  oz*ml.L82,1),2);
dwdylse=fftshift(fftshift(v*ml.L91+  oz*ml.L92,1),2);

dudzlse=fftshift(fftshift(v*ml.L101 +oz*ml.L102,1),2);
dvdzlse=fftshift(fftshift(v*ml.L111 +oz*ml.L112,1),2);
dwdzlse=fftshift(fftshift(v*ml.L121 +oz*ml.L122,1),2);


%%
yp=yCheb+1;
[X,Z,Y]=meshgrid(xp,zp,yp(ny/2+1:end));
fn=sprintf('../data/velgrad_voz_field_lseQ4_j_%03d.mat',jcond);
m=matfile(fn,'Writable',true);
m.ozcond=oz;
m.vcond=v;
m.u=ulse;
m.v=vlse;
m.w=wlse;
m.dudx=dudxlse;
m.dvdx=dvdxlse;
m.dwdx=dwdxlse;
m.dudy=dudylse;
m.dvdy=dvdylse;
m.dwdy=dwdylse;
m.dudz=dudzlse;
m.dvdz=dvdzlse;
m.dwdz=dwdzlse;
m.X=X;
m.Y=Y;
m.Z=Z;
%m.fxQ2=fxQ2;
%m.fxQ4=fxQ4;


