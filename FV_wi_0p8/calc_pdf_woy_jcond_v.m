clear 
close all
Nx=512;
Nz=384;
Ny=220;
jcond=188;
tstart=10000;
tend=108000;
tstep=1000;
nf=(tend-tstart)/tstep+1;
nbin=256;
% wmin=-0.1;
% wmax=0.1;
% oymin=-5;
% oymax=5;
wmin=-0.3;
wmax=0.3;
oymin=-15;
oymax=15;
densityn=zeros(nbin,nbin);
densityp=zeros(nbin,nbin);
%load('lambda_stats.mat')
for time=tstart:tstep:tend
        time
        fvel=sprintf("velfields_%07d.mat",time);
        m=matfile(fvel);
	fvg=sprintf("velgrad_%07d.mat",time);
            mg=matfile(fvg);

	oy=mg.dudz(:,:,jcond)-mg.dwdx(:,:,jcond);
	oy=oy-mean(oy,'all');
	w=m.wfield(:,:,jcond)-mean( m.wfield(:,:,jcond),'all');
	v=m.vfield(:,:,jcond);	
	oyp=oy(v>0);
	oyn=oy(v<0);
	wp=w(v>0);
	wn=w(v<0);
        %oyl=reshape( oy ,Nx*Nz,1 );
        %wl=reshape( m.wfield(:,:,jcond)-mean( m.wfield(:,:,jcond),'all') ,Nx*Nz,1 );
        woyp=[wp,oyp];
	woyn=[wn,oyn];

        [bandwidth,densityinstp,wbin,oybin]=kde2d(woyp,nbin,[wmin oymin],[wmax oymax]);
	[bandwidth,densityinstn,wbin,oybin]=kde2d(woyn,nbin,[wmin oymin],[wmax oymax]);
        densityp=densityp+densityinstp;
	densityn=densityn+densityinstn;
end
densityn=densityn./nf;
densityp=densityp./nf;

fn=sprintf('jointpdf_woy_v_j_%03d.mat',jcond);
mf=matfile(fn,"Writable",true);
mf.densityp=densityp;
mf.densityn=densityn;
mf.bandwidth=bandwidth;
mf.wbin=wbin;
mf.oybin=oybin;
