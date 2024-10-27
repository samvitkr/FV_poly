close all
clear
jcond=194;

nx=512;
nz=384;
ny=220;
ulse=zeros(nz,nx,ny/2);
vlse=zeros(nz,nx,ny/2);
wlse=zeros(nz,nx,ny/2);
fxlse=zeros(nz,nx,ny/2);

lx=4*pi;
lz=2*pi;
xp=lx*[0:nx-1]/nx-lx/2;
zp=lz*[0:nz-1]/nz-lz/2;

load('../data/mean_profiles.mat')

%U=reshape(Um,[1 1 ny]);
%dUdy=reshape(dUdy,[1 1 ny]);
%viscm=reshape(viscm,[ 1 1 ny]);
%U=U(:,:,111:end);
%dUdy=dUdy(:,:,111:end);
%viscm=viscm(:,:,111:end);

fp=sprintf('../data/jointpdf_uv_j_%03d.mat',jcond);
mp=matfile(fp);

fn=sprintf('../data/velgrad_lse_j_%03d.mat',jcond);
ml=matfile(fn);
%%

u=mp.u2;
v=mp.v2;

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

fxlse=fftshift(fftshift(u*ml.L131+v*ml.L132,1),2);

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
m.fx=fxlse;

%%

u=mp.u4;	
v=mp.v4;

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

fxlse=fftshift(fftshift(u*ml.L131+v*ml.L132,1),2);


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
m.fx=fxlse;


