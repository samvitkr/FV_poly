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

fn=sprintf('../data/lsefilter_j_%03d.mat',jcond)
mf=matfile(fn);
load('../data/ygrid.mat')
%U=reshape(Um,[1 1 ny]);
%dUdy=reshape(dUdy,[1 1 ny]);
%viscm=reshape(viscm,[ 1 1 ny]);
%U=U(:,:,111:end);
%dUdy=dUdy(:,:,111:end);
%viscm=viscm(:,:,111:end);



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

fxlse=fftshift(fftshift(u*ml.L131+v*ml.L132,1),2);

yp=yCheb+1;
[X,Z,Y]=meshgrid(xp,zp,yp(ny/2+1:end));
fn=sprintf('../data/velgradfield_dfil_lseQ2_j_%03d.mat',jcond);
m=matfile(fn,'Writable',true);
m.ucond=u;
m.vcond=v;
m.u=ifft2(fft2(ulse).*mf.dfil);
m.v=ifft2(fft2(vlse).*mf.dfil);
m.w=ifft2(fft2(wlse).*mf.dfil);
m.dudx=ifft2(fft2(dudxlse).*mf.dfil);
m.dvdx=ifft2(fft2(dvdxlse).*mf.dfil);
m.dwdx=ifft2(fft2(dwdxlse).*mf.dfil);
m.dudy=ifft2(fft2(dudylse).*mf.dfil);
m.dvdy=ifft2(fft2(dvdylse).*mf.dfil);
m.dwdy=ifft2(fft2(dwdylse).*mf.dfil);
m.dudz=ifft2(fft2(dudzlse).*mf.dfil);
m.dvdz=ifft2(fft2(dvdzlse).*mf.dfil);
m.dwdz=ifft2(fft2(dwdzlse).*mf.dfil);
m.X=X;
m.Y=Y;
m.Z=Z;
m.fx=ifft2(fft2(fxlse).*mf.dfil);


fnu=sprintf('../data/velgradfield_ufil_lseQ2_j_%03d.mat',jcond);
mu=matfile(fnu,'Writable',true);
mu.ucond=u;
mu.vcond=v;
mu.u=ifft2(fft2(ulse).*mf.ufil);
mu.v=ifft2(fft2(vlse).*mf.ufil);
mu.w=ifft2(fft2(wlse).*mf.ufil);
mu.dudx=ifft2(fft2(dudxlse).*mf.ufil);
mu.dvdx=ifft2(fft2(dvdxlse).*mf.ufil);
mu.dwdx=ifft2(fft2(dwdxlse).*mf.ufil);
mu.dudy=ifft2(fft2(dudylse).*mf.ufil);
mu.dvdy=ifft2(fft2(dvdylse).*mf.ufil);
mu.dwdy=ifft2(fft2(dwdylse).*mf.ufil);
mu.dudz=ifft2(fft2(dudzlse).*mf.ufil);
mu.dvdz=ifft2(fft2(dvdzlse).*mf.ufil);
mu.dwdz=ifft2(fft2(dwdzlse).*mf.ufil);
mu.X=X;
mu.Y=Y;
mu.Z=Z;
mu.fx=ifft2(fft2(fxlse).*mf.ufil);

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

fxlse=fftshift(fftshift(u*ml.L131+v*ml.L132,1),2);


%%
yp=yCheb+1;
[X,Z,Y]=meshgrid(xp,zp,yp(ny/2+1:end));
fn=sprintf('../data/velgradfield_dfil_lseQ4_j_%03d.mat',jcond);
m=matfile(fn,'Writable',true);
m.ucond=u;
m.vcond=v;
m.u=ifft2(fft2(ulse).*mf.dfil);
m.v=ifft2(fft2(vlse).*mf.dfil);
m.w=ifft2(fft2(wlse).*mf.dfil);
m.dudx=ifft2(fft2(dudxlse).*mf.dfil);
m.dvdx=ifft2(fft2(dvdxlse).*mf.dfil);
m.dwdx=ifft2(fft2(dwdxlse).*mf.dfil);
m.dudy=ifft2(fft2(dudylse).*mf.dfil);
m.dvdy=ifft2(fft2(dvdylse).*mf.dfil);
m.dwdy=ifft2(fft2(dwdylse).*mf.dfil);
m.dudz=ifft2(fft2(dudzlse).*mf.dfil);
m.dvdz=ifft2(fft2(dvdzlse).*mf.dfil);
m.dwdz=ifft2(fft2(dwdzlse).*mf.dfil);
m.X=X;
m.Y=Y;
m.Z=Z;
m.fx=ifft2(fft2(fxlse).*mf.dfil);


fnu=sprintf('../data/velgradfield_ufil_lseQ4_j_%03d.mat',jcond);
mu=matfile(fnu,'Writable',true);
mu.ucond=u;
mu.vcond=v;
mu.u=ifft2(fft2(ulse).*mf.ufil);
mu.v=ifft2(fft2(vlse).*mf.ufil);
mu.w=ifft2(fft2(wlse).*mf.ufil);
mu.dudx=ifft2(fft2(dudxlse).*mf.ufil);
mu.dvdx=ifft2(fft2(dvdxlse).*mf.ufil);
mu.dwdx=ifft2(fft2(dwdxlse).*mf.ufil);
mu.dudy=ifft2(fft2(dudylse).*mf.ufil);
mu.dvdy=ifft2(fft2(dvdylse).*mf.ufil);
mu.dwdy=ifft2(fft2(dwdylse).*mf.ufil);
mu.dudz=ifft2(fft2(dudzlse).*mf.ufil);
mu.dvdz=ifft2(fft2(dvdzlse).*mf.ufil);
mu.dwdz=ifft2(fft2(dwdzlse).*mf.ufil);
mu.X=X;
mu.Y=Y;
mu.Z=Z;
mu.fx=ifft2(fft2(fxlse).*mf.ufil);
