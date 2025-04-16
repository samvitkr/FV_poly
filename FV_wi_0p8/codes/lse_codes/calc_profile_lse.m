close all
clear
jcond=156;
ltq=22;
lt=0.005;
%fth=0.1;
%ft=sprintf("velfield_lse_voz_j_%d.mat",jcond)
%ft=sprintf("velfield_lse_vwoy_j_%d.mat",jcond)
load('../data/ygrid.mat')
%ft =sprintf('../data/velgradfield_lseQ2_j_%03d.mat',jcond);
ft=sprintf('../data/velgradfx_voz_field_lseQ2_j_156.mat')
%ft=sprintf('../data/velgradfield_lseQ2_voz_j_156.mat')
m=matfile(ft)
%ftu=sprintf('../data/velgradfield_lseQ4_j_%03d.mat',jcond);
ftu=sprintf('../data/velgradfx_voz_field_lseQ4_j_156.mat')
%ftu=sprintf('../data/velgradfield_lseQ4_voz_j_156.mat')

mu=matfile(ftu)

qrmsd=rms(m.Q,[1 2]);
qrmsu=rms(mu.Q,[1 2]);

nl2=(m.v).*(m.dvdx-m.dudz)-(m.w).*(m.dudz-m.dwdx);
nl4=(mu.v).*(mu.dvdx-mu.dudz)-(mu.w).*(mu.dudz-mu.dwdx);
voz2=(m.v).*(m.dvdx-m.dudz);
woy2=(m.w).*(m.dudz-m.dwdx);
voz4=(mu.v).*(mu.dvdx-mu.dudz);
woy4=(mu.w).*(mu.dudz-mu.dwdx);

%nl2rms=(rms(nl2,[1 2]));
%nl4rms=(rms(nl4,[1 2]));

nl2m=squeeze(mean(nl2,[1 2]));
nl4m=squeeze(mean(nl4,[1 2]));
voz2m=squeeze(mean(voz2,[1 2]));
voz4m=squeeze(mean(voz4,[1 2]));
woy2m=squeeze(mean(woy2,[1 2]));
woy4m=squeeze(mean(woy4,[1 2]));

nl2q=squeeze(mean(nl2.*(m.Q>0),[1 2]));
nl4q=squeeze(mean(nl4.*(mu.Q>0),[1 2]));
nl2qrms=squeeze(mean(nl2.*(m.Q./qrmsd>ltq),[1 2]));
nl4qrms=squeeze(mean(nl4.*(mu.Q./qrmsu>ltq),[1 2]));



% nl2d=squeeze(mean(nl2.*(nl2<0),[1 2]));
% nl2u=squeeze(mean(nl2.*(nl2>0),[1 2]));
% nl4d=squeeze(mean(nl4.*(nl4<0),[1 2]));
% nl4u=squeeze(mean(nl4.*(nl4>0),[1 2]));
% nl2r=squeeze(mean(nl2.*(abs(nl2)./nl2rms>1),[1 2]));
% nl4r=squeeze(mean(nl4.*(abs(nl4)./nl4rms>1),[1 2]));

fpu=sprintf('../data/velgradfx_voz_lse_profile_j_%03d.mat',jcond);
mpu=matfile(fpu,'Writable',true)
mpu.nl2m=nl2m;
mpu.nl4m=nl4m;
mpu.voz2m=voz2m;
mpu.voz4m=voz4m;
mpu.woy2m=woy2m;
mpu.woy4m=woy4m;
mpu.nl2q=nl2q; 
mpu.nl4q=nl4q;
mpu.nl2qrms=nl2qrms;                       
mpu.nl4qrms=nl4qrms;
% mpu.nl2d=nl2d;
% mpu.nl2u=nl2u;
% mpu.nl4d=nl4d;
% mpu.nl4u=nl4u;
% mpu.nl2r=nl2r;
% mpu.nl4r=nl4r;

