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
vmin=-0.1;
vmax=0.1;
ozmin=-7;
ozmax=7;
density=zeros(nbin,nbin);

%load('lambda_stats.mat')
for time=tstart:tstep:tend
        time
        fvel=sprintf("velfields_%07d.mat",time);
        m=matfile(fvel);
	fvg=sprintf("velgrad_%07d.mat",time);
            mg=matfile(fvg);

	oz=mg.dvdx(:,:,jcond)-mg.dudy(:,:,jcond);
	oz=oz-mean(oz,'all');
	
        ozl=reshape( oz ,Nx*Nz,1 );
        vl=reshape( m.vfield(:,:,jcond) ,Nx*Nz,1 );
        voz=[vl,ozl];
        [bandwidth,densityinst,vbin,ozbin]=kde2d(voz,nbin,[vmin ozmin],[vmax ozmax]);
        density=density+densityinst;
end
density=density./nf;

fn=sprintf('jointpdf_voz_j_%03d.mat',jcond);
mf=matfile(fn,"Writable",true);
mf.density=density;
mf.bandwidth=bandwidth;
mf.vbin=vbin;
mf.ozbin=ozbin;
