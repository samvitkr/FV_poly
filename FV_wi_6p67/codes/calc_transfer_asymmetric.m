%addpath '/home/skumar67/data-geyink1/skumar67/FV_wi_0p8'
tstart=00000;
tend=60000;
tstep=1000;
Nx=512;
Ny=220;
Nz=384;
load('ygrid.mat')
load('dpdx.mat')
%voz=zeros(Nz,Nx,Ny);
%woy=zeros(Nz,Nx,Ny);
%visc=zeros(Nz,Nx,Ny);
%poly=zeros(Nz,Nx,Ny);
%U=voz;
%vozF=voz;
%woyF=woy;
%uu=voz;
%vv=voz;
%ww=voz;
%uv=voz;
%uw=voz;
%vw=voz;
nt=(tend-tstart)/tstep+1;
syz_0_pi=zeros(nt,1);
syz_pi_2pi=zeros(nt,1);
sydiff=zeros(nt,1);
sydiffp=zeros(nt,1);
ut=ut_ts(1+tstart:tstep:tend+1);
for time=tstart:tstep:tend
	
	ft=sprintf("transferfields_%07d.mat",time)
	mt=matfile(ft);
	
	%voz=voz+mt.voz;
	%woy=woy+mt.woy;
	%vozF=vozF+mt.vozF;
        %woyF=woyF+mt.woyF;
	%visc=visc+mt.visc;
 	%poly=poly+mt.poly;
	syz=mt.voz-mt.woy+mt.visc+mt.poly;
	syz1=squeeze(mean(mean( syz(1:Nz/2,:,:) ,1),2));
	syz2=squeeze(mean(mean( syz(1+Nz/2:end,:,:) ,1),2));
	%size(yCheb)
	%size(syz1)
	%size(syz2)
	%time
	syz_0_pi(time+1)=0.25*trapz( flip(yCheb) , flip(syz1));
	syz_pi_2pi(time+1)=0.25*trapz( flip(yCheb) , flip(syz2));


%        fv=sprintf("velfields_%07d.mat",time)
%        mv=matfile(fv);
%	U=U+mv.ufield;
%	uu=uu+(mv.ufield).*(mv.ufield);
%	vv=vv+(mv.vfield).*(mv.vfield);
%	ww=ww+(mv.wfield).*(mv.wfield);
%	uv=uv+(mv.ufield).*(mv.vfield);
%	uw=uw+(mv.ufield).*(mv.wfield);
%	vw=vw+(mv.vfield).*(mv.wfield);
end
mn=matfile('transferfields_asym_mean.mat','Writable',true)
t=[tstart:tstep:tend]';
mn.t=t;
mn.ut=ut;
mn.syz_0_pi= syz_0_pi(t+1);
mn.syz_pi_2pi= syz_pi_2pi(t+1);
%mn.syzdiff=syz_pi_2pi-syz_0_pi;
%mn.syzdiffp=-(mn.syzdiff)./(ut.^2);
%[valp indp]=max(abs (mn.syzdiffp) );
%[val ind]=max(abs (mn.syzdiff) );
%mn.val=val;
%mn.ind=t(ind);
%mn.valp=valp;
%mn.indp=t(indp);
