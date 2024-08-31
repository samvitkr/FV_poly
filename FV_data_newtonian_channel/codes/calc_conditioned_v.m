Ny=220;
re=4667;
ret=231.8;
ut = ret/re;
dnu=1/ret;
njstart=1;
njend=Ny;
nj=njend-njstart+1;
edges=[-inf,0,inf];
ne=length(edges);
%p=zeros(Ny/2,ne-1);
p=zeros(nj,ne-1);
pvoz=p;
pwoy=p;
ppoly=p;
pv=p;
pvisc=p;

tstart=00000;
tend=40000;
tstep=1000;
nt=(tend-tstart)/tstep+1;
%load('lambda_stats.mat')
for time=tstart:tstep:tend
	time
	fvel=sprintf("velfields_%07d.mat",time);
	mvel=matfile(fvel);
	v=mvel.vfield(:,:,njstart:njend);
	ft=sprintf("transferfields_%07d.mat",time)
	mcv=matfile(ft);
	voz=mcv.voz(:,:,:);%mcv.convective_x(njstart:njend,:,:);
    	woy=mcv.woy(:,:,:);%mcv.viscous_x(njstart:njend,:,:);
	visc=mcv.visc;	
	%%
	
	for i =1:ne-1
		i

		l=v;
                [counts, ~, bins] = histcounts(l(:,:,:), [edges(i) edges(i+1)]);
                size(l)
                size(p)
                size(mean(mean(bins,1),2))
                pv(:,i)  =pv(:,i)  +squeeze(mean(mean(bins,1),2));
                pvoz(:,i)=pvoz(:,i)+squeeze(mean(mean(bins.*voz(:,:,:),1),2));%./p(:,i);
                pwoy(:,i)=pwoy(:,i)+squeeze(mean(mean(bins.*woy(:,:,:),1),2));%./p(:,i);
                pvisc(:,i)=pvisc(:,i)+squeeze(mean(mean(bins.*visc(:,:,:),1),2));%
	end
end
pv=pv./nt;
pvoz=pvoz/nt;
pwoy=pwoy/nt;
pvisc=pvisc/nt;

mae=matfile('conditioned_v','Writable',true);
mae.edges=edges;
mae.pv=pv;
mae.pvoz=pvoz;
mae.pwoy=pwoy;
mae.pvisc=pvisc;
