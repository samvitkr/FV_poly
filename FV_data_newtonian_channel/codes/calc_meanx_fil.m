close all
clear
nx=512;
nz=384;
ny=220;
lx=4*pi;
lz=2*pi;
xp=lx*[0:nx-1]/nx-lx/2;
zp=lz*[0:nz-1]/nz-lz/2;
jcond=194;
%lt=-0.001;
%lt=0.22;vim ca
%lt=-0.005;
%fth=0.1;
%ft=sprintf("velfield_lse_voz_j_%d.mat",jcond)
%ft=sprintf("velfield_lse_vwoy_j_%d.mat",jcond)
load('../data/ygrid.mat')
mm=matfile('../data/mean_profiles.mat')
ut=mm.ret/mm.re;
%dUdy=reshape(mm.dUdy,[1,1,ny]);
%dUdy=dUdy(1,1,111:end);
yp=yCheb(111:end)'+1;
[Z,Y]=meshgrid(zp,yp);

%ft =sprintf('../data/velgradfield_dfil_lseQ4_j_%03d.mat',jcond);
%ft =sprintf('../data/velgradfield_lsevp_j_%03d.mat',jcond);
ft=sprintf('../data/velgradfield_dfil_lseQ2_j_%03d.mat',jcond')
m=matfile(ft,'Writable',true)
%ftu =sprintf('../data/velgradfield_lsevn_j_%03d.mat',jcond);
%ftu=sprintf('../data/velgrad_voz_field_lseQ4_j_156.mat')
ftu=sprintf('../data/velgradfield_ufil_lseQ4_j_%03d.mat',jcond')
mu=matfile(ftu,'Writable',true)

fx=sprintf('../data/velgrad_meanx_lse_fil_Q2dQ4u_j_%03d.mat',jcond);
mx=matfile(fx,'Writable',true)

ltd=m.ltot;
vd=squeeze(mean( (m.v).*(m.lambda2tot<ltd),2))';
wd=squeeze(mean( (m.w).*(m.lambda2tot<ltd),2))';
ozd=squeeze(mean( (m.dvdx-m.dudy-m.dUdym).*(m.lambda2tot<ltd),2))';
oyd=squeeze(mean( (m.dudz-m.dwdx).*(m.lambda2tot<ltd),2))';
vozd=squeeze(mean( (m.v).*(m.dvdx-m.dudy-m.dUdym).*(m.lambda2tot<ltd),2))';
woyd=squeeze(mean( (m.w).*(m.dudz-m.dwdx).*(m.lambda2tot<ltd),2))';
ld=squeeze(mean( (m.lambda2tot).*(m.lambda2tot<ltd),2))';

ltu=mu.ltot;
vu=squeeze(mean( (mu.v).*(mu.lambda2tot<ltu),2))';
wu=squeeze(mean( (mu.w).*(mu.lambda2tot<ltu),2))';
ozu=squeeze(mean( (mu.dvdx-mu.dudy-m.dUdym).*(mu.lambda2tot<ltu),2))';
oyu=squeeze(mean( (mu.dudz-mu.dwdx).*(mu.lambda2tot<ltu),2))';
vozu=squeeze(mean( (mu.v).*(mu.dvdx-mu.dudy-m.dUdym).*(mu.lambda2tot<ltu),2))';
woyu=squeeze(mean( (mu.w).*(mu.dudz-mu.dwdx).*(mu.lambda2tot<ltu),2))';
lu=squeeze(mean( (mu.lambda2tot).*(mu.lambda2tot<ltu),2))';

mx.vd=vd;
mx.vu=vu;
mx.ozd=ozd;
mx.ozu=ozu;
mx.wd=wd;
mx.wu=wu;
mx.oyd=oyd;
mx.oyu=oyu;
mx.vozd=vozd;
mx.vozu=vozu;
mx.woyd=woyd;
mx.woyu=woyu;
mx.ld=ld;
mx.lu=lu;
mx.Z=Z;
mx.Y=Y;

%%
nld=vozd-woyd;
nlu=vozu-woyu;
nldm=mean(nld,2);
nlum=mean(nlu,2);

hold on
plot(-nldm,yp,'-r')
plot(-nlum,yp,'-b')
hold off
grid on



%%
%  pcolor(Z,Y,woyd)
%  shading flat
%  axis equal
% clim([-1e-5 1e-5])
