clear 
close all
Nx=512;
Nz=384;
Ny=220;
jcond=156;
jc=jcond-110;
tstart=0000;
tend=100000;
tstep=1000;
nf=(tend-tstart)/tstep+1;
nbin=256;
umin=-0.25;
umax=0.25;
vmin=-0.25;
vmax=0.25;
densityd=zeros(nbin,nbin);
densityu=zeros(nbin,nbin);
mf=matfile('../data/filter.mat')
dfil=mf.dfil;
ufil=mf.ufil;
%load('lambda_stats.mat')
for time=tstart:tstep:tend
        time
        fvel=sprintf("../data/velfields_%07d.mat",time);
         m=matfile(fvel);
	 uldf= ifft2(fft2(m.ufield(:,:,jcond)).*dfil(:,:,jc),'symmetric');
	 uluf= ifft2(fft2(m.ufield(:,:,jcond)).*ufil(:,:,jc),'symmetric');
	 vldf= ifft2(fft2(m.vfield(:,:,jcond)).*dfil(:,:,jc),'symmetric');
         vluf= ifft2(fft2(m.vfield(:,:,jcond)).*ufil(:,:,jc),'symmetric');
         uld=reshape( uldf-mean(uldf,'all'),Nx*Nz,1 );
         vld=reshape( vldf-mean(vldf,'all'),Nx*Nz,1 );
	 ulu=reshape( uluf-mean(uluf,'all'),Nx*Nz,1 );
         vlu=reshape( vluf-mean(vluf,'all'),Nx*Nz,1 );
        uvd=[uld,vld];
	uvu=[ulu,vlu];
        [bandwidth,densityinst,ubin,vbin]=kde2d(uvd,nbin,[umin vmin],[umax vmax]);
        densityd=densityd+densityinst;
	[bandwidth,densityinst,ubin,vbin]=kde2d(uvu,nbin,[umin vmin],[umax vmax]);
        densityu=densityu+densityinst;
end
densityd=densityd./nf;
densityu=densityu./nf;

fn=sprintf('../data/jointpdf_uvfil_j_%03d.mat',jcond);
mf=matfile(fn,"Writable",true);
mf.densityd=densityd;
mf.densityu=densityu;
mf.bandwidth=bandwidth;
mf.ubin=ubin;
mf.vbin=vbin;
