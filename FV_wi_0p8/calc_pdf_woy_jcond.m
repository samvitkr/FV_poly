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
wmin=-0.1;
wmax=0.1;
oymin=-5;
oymax=5;
density=zeros(nbin,nbin);

%load('lambda_stats.mat')
for time=tstart:tstep:tend
        time
        fvel=sprintf("velfields_%07d.mat",time);
        m=matfile(fvel);
	fvg=sprintf("velgrad_%07d.mat",time);
            mg=matfile(fvg);

	oy=mg.dudz(:,:,jcond)-mg.dwdx(:,:,jcond);
	oy=oy-mean(oy,'all');
	
        oyl=reshape( oy ,Nx*Nz,1 );
        wl=reshape( m.wfield(:,:,jcond)-mean( m.wfield(:,:,jcond),'all') ,Nx*Nz,1 );
        woy=[wl,oyl];
        [bandwidth,densityinst,wbin,oybin]=kde2d(woy,nbin,[wmin oymin],[wmax oymax]);
        density=density+densityinst;
end
density=density./nf;

fn=sprintf('jointpdf_woy_j_%03d.mat',jcond);
mf=matfile(fn,"Writable",true);
mf.density=density;
mf.bandwidth=bandwidth;
mf.wbin=wbin;
mf.oybin=oybin;