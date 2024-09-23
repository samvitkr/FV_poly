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
ft= sprintf('../data/velgradfx_voz_field_lseQ2_j_156.mat')
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
visc2=m.fx;
visc4=mu.fx;
%nl2rms=(rms(nl2,[1 2]));
%nl4rms=(rms(nl4,[1 2]));

nl2m=squeeze(mean(nl2,[1 2]));
nl4m=squeeze(mean(nl4,[1 2]));
voz2m=squeeze(mean(voz2,[1 2]));
voz4m=squeeze(mean(voz4,[1 2]));
woy2m=squeeze(mean(woy2,[1 2]));
woy4m=squeeze(mean(woy4,[1 2]));
visc2m=squeeze( mean( visc2, [1 2] ) );
visc4m=squeeze( mean( visc4, [1 2] ) );


nl2q=squeeze(mean(nl2.*(m.Q>lt),[1 2]));
nl4q=squeeze(mean(nl4.*(mu.Q>lt),[1 2]));
nl2qrms=squeeze(mean(nl2.*(m.Q./qrmsd>lt),[1 2]));
nl4qrms=squeeze(mean(nl4.*(mu.Q./qrmsu>lt),[1 2]));

visc2q=squeeze(mean(visc2.*(m.Q>lt),[1 2]));
visc4q=squeeze(mean(visc4.*(mu.Q>lt),[1 2]));
visc2qrms=squeeze(mean(visc2.*(m.Q./qrmsd>lt),[1 2]));
visc4qrms=squeeze(mean(visc4.*(mu.Q./qrmsu>lt),[1 2]));

voz2q=squeeze(mean(voz2.*(m.Q>lt),[1 2]));
voz4q=squeeze(mean(voz4.*(mu.Q>lt),[1 2]));
woy2q=squeeze(mean(woy2.*(m.Q>lt),[1 2]));
woy4q=squeeze(mean(woy4.*(mu.Q>lt),[1 2]));

voz2qrms=squeeze(mean(voz2.*(m.Q./qrmsd>ltq),[1 2]));
voz4qrms=squeeze(mean(voz4.*(mu.Q./qrmsu>ltq),[1 2]));
woy2qrms=squeeze(mean(woy2.*(m.Q./qrmsd>ltq),[1 2]));
woy4qrms=squeeze(mean(woy4.*(mu.Q./qrmsu>ltq),[1 2]));



% nl2d=squeeze(mean(nl2.*(nl2<0),[1 2]));
% nl2u=squeeze(mean(nl2.*(nl2>0),[1 2]));
% nl4d=squeeze(mean(nl4.*(nl4<0),[1 2]));
% nl4u=squeeze(mean(nl4.*(nl4>0),[1 2]));
% nl2r=squeeze(mean(nl2.*(abs(nl2)./nl2rms>1),[1 2]));
% nl4r=squeeze(mean(nl4.*(abs(nl4)./nl4rms>1),[1 2]));

fpu=sprintf('../data/velgrad_voz_lse_profile_j_%03d.mat',jcond);
%fpu=sprintf('../data/velgrad_lse_profile_j_%03d.mat',jcond);

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

mpu.visc2q=visc2q; 
mpu.visc4q=visc4q;
mpu.visc2qrms=visc2qrms;                       
mpu.visc4qrms=visc4qrms;

mpu.voz2q=voz2q; 
mpu.voz4q=voz4q;
mpu.woy2q=woy2q; 
mpu.woy4q=woy4q;

mpu.voz2qrms=voz2qrms;                       
mpu.voz4qrms=voz4qrms;
mpu.woy2qrms=woy2qrms;                       
mpu.woy4qrms=woy4qrms;
