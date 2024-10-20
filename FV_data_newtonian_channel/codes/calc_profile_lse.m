close all
clear
jcond=194;
ltq=22;
lt=-0.005;
ny=220;
%fth=0.1;
%ft=sprintf("velfield_lse_voz_j_%d.mat",jcond)
%ft=sprintf("velfield_lse_vwoy_j_%d.mat",jcond)
load('../data/ygrid.mat')
%ft =sprintf('../data/velgradfield_dfil_lseQ4_j_%03d.mat',jcond);
%ft= sprintf('../data/velgrad_voz_field_lseQ2_j_156.mat')
ft=sprintf('../data/velgradfield_lsevp_j_%03d.mat',jcond)
m=matfile(ft)
%ftu=sprintf('../data/velgradfield_ufil_lseQ4_j_%03d.mat',jcond);
%ftu=sprintf('../data/velgrad_voz_field_lseQ4ozp_j_156.mat')
ftu=sprintf('../data/velgradfield_lsevn_j_%03d.mat',jcond)
mu=matfile(ftu)
mp=matfile('../data/mean_profiles.mat');
dUdy=reshape(mp.dUdy,[1 1 ny]);
%qrmsd=rms(m.Q,[1 2]);
%qrmsu=rms(mu.Q,[1 2]);

nl2=(m.v).*(m.dvdx-m.dudz)-(m.w).*(m.dudz-m.dwdx);
nl4=(mu.v).*(mu.dvdx-mu.dudz)-(mu.w).*(mu.dudz-mu.dwdx);
voz2=(m.v).*(m.dvdx-m.dudz);
woy2=(m.w).*(m.dudz-m.dwdx);
voz4=(mu.v).*(mu.dvdx-mu.dudz);
woy4=(mu.w).*(mu.dudz-mu.dwdx);
visc2=m.fx;
visc4=mu.fx;
vOz2=(m.v).*(-dUdy(1,1,111:end,1));
vOz4=(mu.v).*(-dUdy(1,1,111:end,1));


nl2m=squeeze(mean(nl2,[1 2]));
nl4m=squeeze(mean(nl4,[1 2]));
voz2m=squeeze(mean(voz2,[1 2]));
voz4m=squeeze(mean(voz4,[1 2]));
woy2m=squeeze(mean(woy2,[1 2]));
woy4m=squeeze(mean(woy4,[1 2]));
visc2m=squeeze( mean( visc2, [1 2] ) );
visc4m=squeeze( mean( visc4, [1 2] ) );
vOz2m=squeeze(mean(vOz2,[1 2]));
vOz4m=squeeze(mean(vOz4,[1 2]));

vol2q=squeeze(mean((m.lambda2tot<lt),[1 2]));
vol4q=squeeze(mean((mu.lambda2tot<lt),[1 2]));
nl2q=squeeze(mean(nl2.*(m.lambda2tot<lt),[1 2]));
nl4q=squeeze(mean(nl4.*(mu.lambda2tot<lt),[1 2]));
visc2q=squeeze(mean(visc2.*(m.lambda2tot<lt),[1 2]));
visc4q=squeeze(mean(visc4.*(mu.lambda2tot<lt),[1 2]));
voz2q=squeeze(mean(voz2.*(m.lambda2tot<lt),[1 2]));
voz4q=squeeze(mean(voz4.*(mu.lambda2tot<lt),[1 2]));
woy2q=squeeze(mean(woy2.*(m.lambda2tot<lt),[1 2]));
woy4q=squeeze(mean(woy4.*(mu.lambda2tot<lt),[1 2]));
vOz2q=squeeze(mean(vOz2.*(m.lambda2tot<lt),[1 2]));
vOz4q=squeeze(mean(vOz4.*(mu.lambda2tot<lt),[1 2]));



vol2qfl=squeeze(mean((m.lambda2<lt),[1 2]));
vol4qfl=squeeze(mean((mu.lambda2<lt),[1 2]));
nl2qfl=squeeze(mean(nl2.*(m.lambda2<lt),[1 2]));
nl4qfl=squeeze(mean(nl4.*(mu.lambda2<lt),[1 2]));
visc2qfl=squeeze(mean(visc2.*(m.lambda2<lt),[1 2]));
visc4qfl=squeeze(mean(visc4.*(mu.lambda2<lt),[1 2]));
voz2qfl=squeeze(mean(voz2.*(m.lambda2<lt),[1 2]));
voz4qfl=squeeze(mean(voz4.*(mu.lambda2<lt),[1 2]));
woy2qfl=squeeze(mean(woy2.*(m.lambda2<lt),[1 2]));
woy4qfl=squeeze(mean(woy4.*(mu.lambda2<lt),[1 2]));
vOz2qfl=squeeze(mean(vOz2.*(m.lambda2<lt),[1 2]));
vOz4qfl=squeeze(mean(vOz4.*(mu.lambda2<lt),[1 2]));


fpu=sprintf('../data/velgrad_v_profile_j_%03d.mat',jcond);

mpu=matfile(fpu,'Writable',true)
mpu.nl2m=nl2m;
mpu.nl4m=nl4m;
mpu.voz2m=voz2m;
mpu.voz4m=voz4m;
mpu.woy2m=woy2m;
mpu.woy4m=woy4m;

mpu.vOz2m=vOz2m;
mpu.vOz4m=vOz4m;

mpu.nl2q=nl2q; 
mpu.nl4q=nl4q;
mpu.visc2q=visc2q; 
mpu.visc4q=visc4q;
mpu.voz2q=voz2q; 
mpu.voz4q=voz4q;
mpu.woy2q=woy2q; 
mpu.woy4q=woy4q;
mpu.vOz2q=vOz2q;
mpu.vOz4q=vOz4q;
mpu.vol2q=vol2q;
mpu.vol4q=vol4q;


mpu.nl2qfl=nl2qfl; 
mpu.nl4qfl=nl4qfl;
mpu.visc2qfl=visc2qfl; 
mpu.visc4qfl=visc4qfl;
mpu.voz2qfl=voz2qfl; 
mpu.voz4qfl=voz4qfl;
mpu.woy2qfl=woy2qfl; 
mpu.woy4qfl=woy4qfl;
mpu.vOz2qfl=vOz2qfl;
mpu.vOz4qfl=vOz4qfl;
mpu.vol2qfl=vol2qfl;
mpu.vol4qfl=vol4qfl;
