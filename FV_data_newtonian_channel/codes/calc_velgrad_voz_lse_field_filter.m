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

lx=4*pi;
lz=2*pi;
xp=lx*[0:nx-1]/nx-lx/2;
zp=lz*[0:nz-1]/nz-lz/2;

fp=sprintf('../data/jointpdf_uv_j_%03d.mat',jcond);
mp=matfile(fp);
load('../data/mean_profiles.mat')
fn=sprintf('../data/filter.mat',jcond)
mf=matfile(fn);
load('../data/ygrid.mat')
U=reshape(Um,[1 1 ny]);
dUdy=reshape(dUdy,[1 1 ny]);
viscm=reshape(viscm,[ 1 1 ny]);
U=U(:,:,111:end);
dUdy=dUdy(:,:,111:end);
viscm=viscm(:,:,111:end);



fn=sprintf('../data/velgrad_lse_j_%03d.mat',jcond);
ml=matfile(fn);
%%
u=mp.u2;
v=mp.v2;

   ulse=fftshift(fftshift(v*ml.L11+  oz*ml.L12,1),2)+U;...+w*ml.L13;
   vlse=fftshift(fftshift(v*ml.L21+  oz*ml.L22,1),2);...+w*ml.L23;
   wlse=fftshift(fftshift(v*ml.L31+  oz*ml.L32,1),2);...+w*ml.L33;

dudxlse=fftshift(fftshift(v*ml.L41+  oz*ml.L42,1),2);
dvdxlse=fftshift(fftshift(v*ml.L51+  oz*ml.L52,1),2);
dwdxlse=fftshift(fftshift(v*ml.L61+  oz*ml.L62,1),2);

dudylse=fftshift(fftshift(v*ml.L71+  oz*ml.L72,1),2)+dUdy;
dvdylse=fftshift(fftshift(v*ml.L81+  oz*ml.L82,1),2);
dwdylse=fftshift(fftshift(v*ml.L91+  oz*ml.L92,1),2);

dudzlse=fftshift(fftshift(v*ml.L101+ oz*ml.L102,1),2);
dvdzlse=fftshift(fftshift(v*ml.L111+ oz*ml.L112,1),2);
dwdzlse=fftshift(fftshift(v*ml.L121+ oz*ml.L122,1),2);

fxlse=fftshift(fftshift(  v*ml.L131+ oz*ml.L132,1),2)+viscm;

yp=yCheb+1;
[X,Z,Y]=meshgrid(xp,zp,yp(ny/2+1:end));
fn=sprintf('../data/velgradfield_dfil_lseQ2_j_%03d.mat',jcond);
m=matfile(fn,'Writable',true);
m.ucond=u;
m.vcond=v;
m.u=ifft2(fft2(ulse).*mf.dfil,'symmetric');
m.v=ifft2(fft2(vlse).*mf.dfil,'symmetric');
m.w=ifft2(fft2(wlse).*mf.dfil,'symmetric');
m.dudx=ifft2(fft2(dudxlse).*mf.dfil,'symmetric');
m.dvdx=ifft2(fft2(dvdxlse).*mf.dfil,'symmetric');
m.dwdx=ifft2(fft2(dwdxlse).*mf.dfil,'symmetric');
m.dudy=ifft2(fft2(dudylse).*mf.dfil,'symmetric');
m.dvdy=ifft2(fft2(dvdylse).*mf.dfil,'symmetric');
m.dwdy=ifft2(fft2(dwdylse).*mf.dfil,'symmetric');
m.dudz=ifft2(fft2(dudzlse).*mf.dfil,'symmetric');
m.dvdz=ifft2(fft2(dvdzlse).*mf.dfil,'symmetric');
m.dwdz=ifft2(fft2(dwdzlse).*mf.dfil,'symmetric');
m.X=X;
m.Y=Y;
m.Z=Z;
m.fx=ifft2(fft2(fxlse).*mf.dfil,'symmetric');


fnu=sprintf('../data/velgradfield_ufil_lseQ2_j_%03d.mat',jcond);
mu=matfile(fnu,'Writable',true);
mu.ucond=u;
mu.vcond=v;
mu.u=ifft2(fft2(ulse).*mf.ufil,'symmetric');
mu.v=ifft2(fft2(vlse).*mf.ufil,'symmetric');
mu.w=ifft2(fft2(wlse).*mf.ufil,'symmetric');
mu.dudx=ifft2(fft2(dudxlse).*mf.ufil,'symmetric');
mu.dvdx=ifft2(fft2(dvdxlse).*mf.ufil,'symmetric');
mu.dwdx=ifft2(fft2(dwdxlse).*mf.ufil,'symmetric');
mu.dudy=ifft2(fft2(dudylse).*mf.ufil,'symmetric');
mu.dvdy=ifft2(fft2(dvdylse).*mf.ufil,'symmetric');
mu.dwdy=ifft2(fft2(dwdylse).*mf.ufil,'symmetric');
mu.dudz=ifft2(fft2(dudzlse).*mf.ufil,'symmetric');
mu.dvdz=ifft2(fft2(dvdzlse).*mf.ufil,'symmetric');
mu.dwdz=ifft2(fft2(dwdzlse).*mf.ufil,'symmetric');
mu.X=X;
mu.Y=Y;
mu.Z=Z;
mu.fx=ifft2(fft2(fxlse).*mf.ufil,'symmetric');

%%
u=mp.u4;	
v=mp.v4;
   ulse=fftshift(fftshift(v*ml.L11+  oz*ml.L12,1),2)+U;...+w*ml.L13;
   vlse=fftshift(fftshift(v*ml.L21+  oz*ml.L22,1),2);...+w*ml.L23;
   wlse=fftshift(fftshift(v*ml.L31+  oz*ml.L32,1),2);...+w*ml.L33;

dudxlse=fftshift(fftshift(v*ml.L41+  oz*ml.L42,1),2);
dvdxlse=fftshift(fftshift(v*ml.L51+  oz*ml.L52,1),2);
dwdxlse=fftshift(fftshift(v*ml.L61+  oz*ml.L62,1),2);

dudylse=fftshift(fftshift(v*ml.L71+  oz*ml.L72,1),2)+dUdy;
dvdylse=fftshift(fftshift(v*ml.L81+  oz*ml.L82,1),2);
dwdylse=fftshift(fftshift(v*ml.L91+  oz*ml.L92,1),2);

dudzlse=fftshift(fftshift(v*ml.L101+ oz*ml.L102,1),2);
dvdzlse=fftshift(fftshift(v*ml.L111+ oz*ml.L112,1),2);
dwdzlse=fftshift(fftshift(v*ml.L121+ oz*ml.L122,1),2);

fxlse=fftshift(fftshift(  v*ml.L131+ oz*ml.L132,1),2)+viscm;

%%
yp=yCheb+1;
[X,Z,Y]=meshgrid(xp,zp,yp(ny/2+1:end));
fn=sprintf('../data/velgradfield_dfil_lseQ4_j_%03d.mat',jcond);
m=matfile(fn,'Writable',true);
m.ucond=u;
m.vcond=v;
m.u=ifft2(fft2(ulse).*mf.dfil,'symmetric');
m.v=ifft2(fft2(vlse).*mf.dfil,'symmetric');
m.w=ifft2(fft2(wlse).*mf.dfil,'symmetric');
m.dudx=ifft2(fft2(dudxlse).*mf.dfil,'symmetric');
m.dvdx=ifft2(fft2(dvdxlse).*mf.dfil,'symmetric');
m.dwdx=ifft2(fft2(dwdxlse).*mf.dfil,'symmetric');
m.dudy=ifft2(fft2(dudylse).*mf.dfil,'symmetric');
m.dvdy=ifft2(fft2(dvdylse).*mf.dfil,'symmetric');
m.dwdy=ifft2(fft2(dwdylse).*mf.dfil,'symmetric');
m.dudz=ifft2(fft2(dudzlse).*mf.dfil,'symmetric');
m.dvdz=ifft2(fft2(dvdzlse).*mf.dfil,'symmetric');
m.dwdz=ifft2(fft2(dwdzlse).*mf.dfil,'symmetric');
m.X=X;
m.Y=Y;
m.Z=Z;
m.fx=ifft2(fft2(fxlse).*mf.dfil,'symmetric');


fnu=sprintf('../data/velgradfield_ufil_lseQ4_j_%03d.mat',jcond);
mu=matfile(fnu,'Writable',true);
mu.ucond=u;
mu.vcond=v;
mu.u=ifft2(fft2(ulse).*mf.ufil,'symmetric');
mu.v=ifft2(fft2(vlse).*mf.ufil,'symmetric');
mu.w=ifft2(fft2(wlse).*mf.ufil,'symmetric');
mu.dudx=ifft2(fft2(dudxlse).*mf.ufil,'symmetric');
mu.dvdx=ifft2(fft2(dvdxlse).*mf.ufil,'symmetric');
mu.dwdx=ifft2(fft2(dwdxlse).*mf.ufil,'symmetric');
mu.dudy=ifft2(fft2(dudylse).*mf.ufil,'symmetric');
mu.dvdy=ifft2(fft2(dvdylse).*mf.ufil,'symmetric');
mu.dwdy=ifft2(fft2(dwdylse).*mf.ufil,'symmetric');
mu.dudz=ifft2(fft2(dudzlse).*mf.ufil,'symmetric');
mu.dvdz=ifft2(fft2(dvdzlse).*mf.ufil,'symmetric');
mu.dwdz=ifft2(fft2(dwdzlse).*mf.ufil,'symmetric');
mu.X=X;
mu.Y=Y;
mu.Z=Z;
mu.fx=ifft2(fft2(fxlse).*mf.ufil,'symmetric');