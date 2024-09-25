clear
close all
jcond=184;
fn=sprintf('../data/jointpdf_uv_j_%03d.mat',jcond);
mf=matfile(fn,"Writable",true);
load(fn)
[mxv,idmax]=min(vbin.*ubin.*density.*(vbin>0),[],'all');
[imax,jmax]=ind2sub(size(density),idmax);
mf.u2=ubin(imax,jmax);
mf.v2=vbin(imax,jmax);

[mxv,idmax]=min(vbin.*ubin.*density.*(vbin<0),[],'all');
[imax,jmax]=ind2sub(size(density),idmax);
mf.u4=ubin(imax,jmax);
mf.v4=vbin(imax,jmax);

