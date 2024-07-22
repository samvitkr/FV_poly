%function[]=print_figure( outfile, file_format )  
%clear


close all
Nx=513;
Ny=220;
Nz=385;
re=4667;
ret=232;
Lx=  4*pi;
Lz = 2*pi;
xp = [0:Nx-1]*Lx/(Nx);
zp=  [0:1:Nz-1]*Lz/(Nz);
load('ygrid.mat');
y = yCheb;
ut = 0.0499;
dnu=1.0006e-3;
nx=2*82;
%ny=130;
nz=1*165;
time=70000;

fvel=sprintf("vel_%03d.mat",time);
fl=sprintf("lambda_f_2Dinertial_%03d",time);
ft=sprintf("Transfer_f_2Dinertial_%03d.mat",time);
ml=matfile(fl)
mcv=matfile(ft);

nystart=ml.jstart;
nyend=ml.jend;
ny=nyend-nystart+1;
[x,z,y] = meshgrid(zp(1:nz),xp(1:nx),yp(nystart:nyend)+1);

mdata=matfile('isosurface_data_lrmsf2D_oz.mat','Writable',true)

lm  =mean(mean( ml.lambda2,3),2);
lrms=rms(rms(ml.lambda2-lm,3),2);
l   =ml.lambda2(:,1:nx,1:nz);

mdata.x=x;
mdata.y=y;
mdata.z=z;
mdata.omegaz=omegaz;
mdata.omegay=omegay;
mdata.omegax=omegax;
mdata.l=l;
mdata.lrms=lrms;
%%
