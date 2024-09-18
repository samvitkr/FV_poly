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
fn=sprintf('../data/velgrad_lse_j_%03d.mat',jcond);
ml=matfile(fn);
%%
%jcond=156
u=-0.1141;
v=0.0671;

   ulse=fftshift(fftshift(u*ml.L11+  v*ml.L12,1),2);...+w*ml.L13;
   vlse=fftshift(fftshift(u*ml.L21+  v*ml.L22,1),2);...+w*ml.L23;
   wlse=fftshift(fftshift(u*ml.L31+  v*ml.L32,1),2);...+w*ml.L33;

dudxlse=fftshift(fftshift(u*ml.L41+  v*ml.L42,1),2);
dvdxlse=fftshift(fftshift(u*ml.L51+  v*ml.L52,1),2);
dwdxlse=fftshift(fftshift(u*ml.L61+  v*ml.L62,1),2);

dudylse=fftshift(fftshift(u*ml.L71+  v*ml.L72,1),2);
dvdylse=fftshift(fftshift(u*ml.L81+  v*ml.L82,1),2);
dwdylse=fftshift(fftshift(u*ml.L91+  v*ml.L92,1),2);

dudzlse=fftshift(fftshift(u*ml.L101+v*ml.L102,1),2);
dvdzlse=fftshift(fftshift(u*ml.L111+v*ml.L112,1),2);
dwdzlse=fftshift(fftshift(u*ml.L121+v*ml.L122,1),2);
%fxlse=u*ml.L41+v*ml.L42+w*ml.L43;

%ulse=(fftshift(fftshift(ulse,1),2));
%vlse=(fftshift(fftshift(vlse,1),2));
%wlse=(fftshift(fftshift(wlse,1),2));

%fxlse=(fftshift(fftshift(fxlse,1),2));
%uQ2(:,:,ny/2+1:end)=ulse;
%vQ2(:,:,ny/2+1:end)=vlse;
%wQ2(:,:,ny/2+1:end)=wlse;
%fxQ2(:,:,ny/2+1:end)=fxlse;

%uQ2(:,:,1:ny/2)=flip(ulse,3);
%vQ2(:,:,1:ny/2)=flip(-vlse,3);
%wQ2(:,:,1:ny/2)=flip(wlse,3);
%fxQ2(:,:,1:ny/2)=flip(fxlse,3);
yp=yCheb+1;
[X,Z,Y]=meshgrid(xp,zp,yp(ny/2+1:end));
fn=sprintf('../data/velgradfield_lseQ2_j_%03d.mat',jcond);
m=matfile(fn,'Writable',true);
m.ucond=u;
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
jcond=156;
u=0.0976;	
v=-0.0576;

   ulse=fftshift(fftshift(u*ml.L11+  v*ml.L12,1),2);...+w*ml.L13;
   vlse=fftshift(fftshift(u*ml.L21+  v*ml.L22,1),2);...+w*ml.L23;
   wlse=fftshift(fftshift(u*ml.L31+  v*ml.L32,1),2);...+w*ml.L33;

dudxlse=fftshift(fftshift(u*ml.L41+  v*ml.L42,1),2);
dvdxlse=fftshift(fftshift(u*ml.L51+  v*ml.L52,1),2);
dwdxlse=fftshift(fftshift(u*ml.L61+  v*ml.L62,1),2);

dudylse=fftshift(fftshift(u*ml.L71+  v*ml.L72,1),2);
dvdylse=fftshift(fftshift(u*ml.L81+  v*ml.L82,1),2);
dwdylse=fftshift(fftshift(u*ml.L91+  v*ml.L92,1),2);

dudzlse=fftshift(fftshift(u*ml.L101+v*ml.L102,1),2);
dvdzlse=fftshift(fftshift(u*ml.L111+v*ml.L112,1),2);
dwdzlse=fftshift(fftshift(u*ml.L121+v*ml.L122,1),2);


%ulse=u*ml.L11+v*ml.L12+w*ml.L13;
%vlse=u*ml.L21+v*ml.L22+w*ml.L23;
%wlse=u*ml.L31+v*ml.L32+w*ml.L33;
%fxlse=u*ml.L41+v*ml.L42+w*ml.L43;
%ulse=(fftshift(fftshift(ulse,1),2));
%vlse=(fftshift(fftshift(vlse,1),2));
%wlse=(fftshift(fftshift(wlse,1),2));
%fxlse=(fftshift(fftshift(fxlse,1),2));
%uQ4(:,:,ny/2+1:end)=ulse;
%vQ4(:,:,ny/2+1:end)=vlse;
%wQ4(:,:,ny/2+1:end)=wlse;

%fxQ4(:,:,ny/2+1:end)=fxlse;
%uQ4(:,:,1:ny/2)=flip(ulse,3);
%vQ4(:,:,1:ny/2)=flip(-vlse,3);
%wQ4(:,:,1:ny/2)=flip(wlse,3);
%fxQ4(:,:,1:ny/2)=flip(fxlse,3);


%%
yp=yCheb+1;
[X,Z,Y]=meshgrid(xp,zp,yp(ny/2+1:end));
fn=sprintf('../data/velgradfield_lseQ4_j_%03d.mat',jcond);
m=matfile(fn,'Writable',true);
m.ucond=u;
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


