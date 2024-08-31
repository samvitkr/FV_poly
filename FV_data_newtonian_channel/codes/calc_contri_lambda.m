Ny=220;
re=4667;
%ret=231.8;
%ut = ret/re;
load('dpdx.mat')
M=matfile('mean_profiles.mat');
vrms=sqrt(M.vv);
ut=mean(ut_ts);
ret=re*ut;
dnu=1/ret;
njstart=1;
njend=Ny;
nj=njend-njstart+1;
edges=[-inf,-1,0,1,inf];
%edges=[-inf,-1,-0.1,0,0.1,1,inf];
ne=length(edges);
%p=zeros(Ny/2,ne-1);
p=zeros(nj,ne-1);
pvoz=p;
pwoy=p;
ppoly=p;
pv=p;
pvisc=p;
poz=p;
tstart=40000;
tend=40000;
tstep=1000;
nt=(tend-tstart)/tstep+1;
%load('lambda_stats.mat')
%for time=tstart:tstep:tend
	%time
	%fvel=sprintf("velfields_%07d.mat",time);
	%mvel=matfile(fvel);
	%v=mvel.vfield(:,:,njstart:njend);
%	v=permute(v,[3 2 1]);
%	v=v./vrms(njstart:njend);
%	v=permute(v,[3 2 1]);

	mcv=matfile('velgrad_transfer_flp_0040000.mat')
	m=matfile('lambdaflp_0040000.mat');
	ml=mean(m.lambda2,[1,2]);
%l=m.lambda2;
	lrms=rms(m.lambda2-ml,[1,2]);
	l=m.lambda2./lrms;
	%ft=sprintf("transferfields_%07d.mat",time)
	%mcv=matfile(ft);
	%fo=sprintf("vortfields_%07d.mat",time)
        %mo=matfile(fo);
	voz=mcv.voz(:,:,:);%mcv.convective_x(njstart:njend,:,:);
    	woy=mcv.woy(:,:,:);%mcv.viscous_x(njstart:njend,:,:);
%	poly=mcv.poly;
%	visc=mcv.visc;
	%oz=mo.omegaZ;	
	%%
	
	for i =1:ne-1
		i

	%	l=v;
                [counts, ~, bins] = histcounts(l(:,:,:), [edges(i) edges(i+1)]);
                size(l)
                size(p)
                size(mean(mean(bins,1),2))
                pv(:,i)  =pv(:,i)  +squeeze(mean(mean(bins,1),2));
                pvoz(:,i)=pvoz(:,i)+squeeze(mean(mean(bins.*voz(:,:,:),1),2));%./p(:,i);
                pwoy(:,i)=pwoy(:,i)+squeeze(mean(mean(bins.*woy(:,:,:),1),2));%./p(:,i);
%		ppoly(:,i)=ppoly(:,i)+squeeze(mean(mean(bins.*poly(:,:,:),1),2));%
%		pvisc(:,i)=pvisc(:,i)+squeeze(mean(mean(bins.*visc(:,:,:),1),2));
%		poz(:,i)=oz:,i)+squeeze(mean(mean(bins.*oz(:,:,:),1),2));
	end
%end
pv=pv./nt;
pvoz=pvoz/nt;
pwoy=pwoy/nt;
%ppoly=ppoly/nt;
%pvisc=pvisc/nt;
mae=matfile('lambdalp_contri','Writable',true);
mae.edges=edges;
mae.pv=pv;
mae.pvoz=pvoz;
mae.pwoy=pwoy;
%mae.ppoly=ppoly;
%mae.pvisc=pvisc;
