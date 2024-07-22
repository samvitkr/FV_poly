%addpath '/home/skumar67/data-geyink1/skumar67/FV_wi_0p8'
tstart=41000;
tend=130000;
tstep=1000;
Nx=512;
Ny=220;
Nz=384;
voz=zeros(Nz,Nx,Ny);
woy=zeros(Nz,Nx,Ny);
visc=zeros(Nz,Nx,Ny);
poly=zeros(Nz,Nx,Ny);
U=voz;
vozF=voz;
woyF=woy;
uu=voz;
vv=voz;
ww=voz;
uv=voz;
uw=voz;
vw=voz;
nt=(tend-tstart)/tstep+1;
for time=tstart:tstep:tend
	
	ft=sprintf("transferfields_%07d.mat",time)
	mt=matfile(ft);
	voz=voz+mt.voz;
	woy=woy+mt.woy;
	vozF=vozF+mt.vozF;
        woyF=woyF+mt.woyF;
	visc=visc+mt.visc;
 	poly=poly+mt.poly;

       fv=sprintf("velfields_%07d.mat",time)
       mv=matfile(fv);
	U=U+mv.ufield;
	uu=uu+(mv.ufield).*(mv.ufield);
	vv=vv+(mv.vfield).*(mv.vfield);
	ww=ww+(mv.wfield).*(mv.wfield);
	uv=uv+(mv.ufield).*(mv.vfield);
	uw=uw+(mv.ufield).*(mv.wfield);
	vw=vw+(mv.vfield).*(mv.wfield);
end
mn=matfile('transferfields_mean.mat','Writable',true)
mn.voz=voz./nt;
mn.woy=woy./nt;
mn.vozF=vozF./nt;
mn.woyF=woyF./nt;
mn.visc=visc./nt;
mn.poly=poly./nt;
mn.U=U./nt;
mn.uu=uu./nt;
mn.vv=vv./nt;
mn.ww=ww./nt;
mn.uv=uv./nt;
mn.uw=uw./nt;
mn.vw=vw./nt;