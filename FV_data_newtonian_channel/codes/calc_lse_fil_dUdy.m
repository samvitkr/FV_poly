Nx=512;
Ny=220;
Nz=384;
jcond=184;

fvg2=sprintf("../data/velgradfield_lseQ2_j_%03d.mat",jcond);
fvg4=sprintf("../data/velgradfield_lseQ4_j_%03d.mat",jcond);
fvgd2=sprintf("../data/velgradfield_dfil_lseQ2_j_%03d.mat",jcond);
fvgu2=sprintf("../data/velgradfield_ufil_lseQ2_j_%03d.mat",jcond);
fvgd4=sprintf("../data/velgradfield_dfil_lseQ4_j_%03d.mat",jcond);
fvgu4=sprintf("../data/velgradfield_ufil_lseQ4_j_%03d.mat",jcond);
fvgn=[fvg2;fvg4;fvgd2; fvgu2; fvgd4; fvgu4];

mf=matfile('../data/filter.mat');
dfilm=mf.dfil(1,1,1:end);
ufilm=mf.ufil(1,1,1:end);
film=[dfilm.*0+1 ;dfilm.*0+1 ;dfilm ; ufilm; dfilm; ufilm;];
mm=matfile('../data/mean_profiles.mat')
dUdy=reshape(mm.dUdy,[1 1 Ny]);
dUdy=dUdy(1,1,111:end);
for nn=1:6
%fvg=fvgoz;
fvg=fvgn(nn);	
fil=film(nn,1,:);
mvg=matfile(fvg,'Writable',true)
mvg.dUdym=dUdy.*fil;
end
