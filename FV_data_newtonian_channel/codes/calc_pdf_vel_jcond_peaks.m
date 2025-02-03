clear 
close all

%jcset=[116,135,187,198,205];
%for jc=1:5
	%clear
	%close all
%	jcond=jcset(jc)
	jcond=171;
	fn=sprintf('../data/jointpdf_uv_j_%03d.mat',jcond);
	mf=matfile(fn,"Writable",true)
	load(fn)
	
	
	%[mxv,idmax]=min(vbin.*ubin.*density.*(vbin>0),[],'all');
	%[imax,jmax]=ind2sub(size(density),idmax);
	%mf.u2=ubin(imax,jmax);
	%mf.v2=vbin(imax,jmax);
	
	%[mxv,idmax]=min(vbin.*ubin.*density.*(vbin<0),[],'all');
	%[imax,jmax]=ind2sub(size(density),idmax);
	%mf.u4=ubin(imax,jmax);
	%mf.v4=vbin(imax,jmax);
	
	 vpdf=sum(vbin.*density,2);
	 [minval,idminval]=min(vpdf);
	 [maxval,idmaxval]=max(vpdf);
	 vminval=vbin(idminval,1)
	 vmaxval=vbin(idmaxval,1)
	 mf.vmin=vminval;
	 mf.vmax=vmaxval;
 %end
