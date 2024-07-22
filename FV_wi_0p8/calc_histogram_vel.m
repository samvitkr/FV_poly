Nx=512;
Nz=384;
jcond=188;
tstart=10000;
tend=108000;
tstep=1000;
nf=(tend-tstart)/tstep+1;
be=0.2*[-1:0.01:1];
sz=length(be)-1;
nval=zeros(1,sz);
%load('lambda_stats.mat')
for time=tstart:tstep:tend
        time
        fvel=sprintf("velfields_%07d.mat",time);
        m=matfile(fvel);
	v=m.vfield(:,:,jcond);
	vl=reshape(v,[1,Nx*Nz]);
	h=histogram( vl,be );
	nval=nval+h.Values;
end
nval=nval./nf;
fn=sprintf('vel_pdf_j_%03d.mat',jcond);
mf=matfile(fn,"Writable",true);
mf.v_BinEdges=be;
mf.v_Values=nval;
