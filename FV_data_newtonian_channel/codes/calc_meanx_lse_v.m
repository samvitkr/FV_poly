jcond=135;
load('../data/ygrid.mat')
fnf=sprintf('../data/velgrad_v_lse_j_%03d.mat',jcond);
mf=matfile(fnf)
nx=512;
nz=384;
ny=220;
lx=4*pi;
lz=2*pi;
xp=lx*[0:nx-1]/nx-lx/2;
zp=lz*[0:nz-1]/nz-lz/2;
mm=matfile('../data/mean_profiles.mat')
ut=mm.ret/mm.re;
dUdy=reshape(mm.dUdy,[1,1,ny]);
dUdy=dUdy(1,1,111:end);
yp=yCheb(111:end)'+1;
[Z,Y]=meshgrid(zp,yp);