close all
clear
jcond=156;
lt=0.007;
%fth=0.1;
%ft=sprintf("velfield_lse_voz_j_%d.mat",jcond)
%ft=sprintf("velfield_lse_vwoy_j_%d.mat",jcond)
load('../data/ygrid.mat')
ft =sprintf('../data/velgradfield_lse_ddfilter_Q2_j_%03d.mat',jcond);
m=matfile(ft)
ftu=sprintf('../data/velgradfield_lse_uufilter_Q2_j_%03d.mat',jcond);
mu=matfile(ftu)

qrmsd=rms(m.Q,[1 2]);
qrmsu=rms(mu.Q,[1 2]);

nld=(m.v).*(m.dvdx-m.dudz)-(m.w).*(m.dudz-m.dwdx);
nlu=(mu.v).*(mu.dvdx-mu.dudz)-(mu.w).*(mu.dudz-mu.dwdx);
nldrms=(rms(nld,[1 2]));
nlurms=(rms(nlu,[1 2]));
nldm=squeeze(mean(nld,[1 2]));
nlum=squeeze(mean(nlu,[1 2]));
nldq=squeeze(mean(nld.*(m.Q>0),[1 2]));
nluq=squeeze(mean(nlu.*(mu.Q>0),[1 2]));
nldqrms=squeeze(mean(nld.*(m.Q./qrmsd>1),[1 2]));
nluqrms=squeeze(mean(nlu.*(mu.Q./qrmsu>1),[1 2]));
nldd=squeeze(mean(nld.*(nld<0),[1 2]));
nldu=squeeze(mean(nld.*(nld>0),[1 2]));
nlud=squeeze(mean(nlu.*(nlu<0),[1 2]));
nluu=squeeze(mean(nlu.*(nlu>0),[1 2]));
nldr=squeeze(mean(nld.*(abs(nld)./nldrms>1),[1 2]));
nlur=squeeze(mean(nlu.*(abs(nlu)./nlurms>1),[1 2]));

fpu=sprintf('../data/velgrad_lse_profile_Q2_j_%03d.mat',jcond);
mpu=matfile(fpu,'Writable',true)
mpu.nldm=nldm;
mpu.nlum=nlum;
mpu.nldq=nldq; 
mpu.nluq=nluq;
mpu.nldqrms=nldqrms;                       
mpu.nluqrms=nluqrms;
mpu.nldd=nldd;
mpu.nldu=nldu;
mpu.nlud=nlud;
mpu.nluu=nluu;
mpu.nldr=nldr;
mpu.nlur=nlur;

