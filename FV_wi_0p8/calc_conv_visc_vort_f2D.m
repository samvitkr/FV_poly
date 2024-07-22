Nx=512;
Ny=220;
Nz=384;
nproc=6;
nzproc=Nz/nproc;
re=4667;
nu=1.0/re;
tstart=50000;
tend=50000;
Nj=Ny;
nt=tend-tstart+1;
lp_lrms=(zeros(Nj,1));
hp_lrms=(zeros(Nj,1));
lp_0=(zeros(Nj,1));
hp_0=(zeros(Nj,1));

volfrac_lp_lrms=(zeros(Nj,1));
volfrac_hp_lrms=(zeros(Nj,1));
volfrac_lp_0=(zeros(Nj,1));
volfrac_hp_0=(zeros(Nj,1));


m=matfile('Transfer_lambda_f2D_0050000.mat','Writable',true)
%        lrmshp=m.lrmshp;
%        lrmslp=m.lrmslp;
for time =tstart:tend
        
       	flhp=sprintf("lambdafhp_%07d",time)
       	mlhp=matfile(flhp);
	lhp=mlhp.lambda2;%./lrmshp;

       	fllp=sprintf("lambdaflp_%07d",time)
       	mllp=matfile(fllp);
	llp=mllp.lambda2;%./lrmslp;

	flt=sprintf("velgrad_transfer_flp_%07d",time)
	mtlp=matfile(flt);
	lp_0=  squeeze(mean(squeeze(mean( ( mtlp.voz-mtlp.woy  ).* ( llp<0 ) )),1));
%	lp_lrms= lp_lrms + mean(squeeze(mean( ( mtlp.voz-mtlp.woy  ).* ( llp<-1 )   ,3 )),2);
	
	volfrac_lp_0= squeeze(mean(squeeze(mean(  ( llp<0 ) )),1));
	%volfrac_lp_lrms=volfrac_lp_lrms + mean(squeeze(mean(  ( llp<-1 )   ,3 )),2);



	fht=sprintf("velgrad_transfer_fhp_%07d",time)
        mthp=matfile(fht);
	hp_0= squeeze(mean(squeeze(mean( ( mthp.voz-mthp.woy  ).* ( lhp<0 ))) ,1)) ;
%        hp_lrms= hp_lrms + mean(squeeze(mean( ( mthp.voz-mthp.woy  ).* ( lhp<-1 ),3 )),2);

	volfrac_hp_0= squeeze(mean(squeeze(mean( ( lhp<0 )   )),1));
%        volfrac_hp_lrms= volfrac_hp_lrms + mean(squeeze(mean( ( lhp<-1 ),3 )),2);

end

m.lp_0=lp_0;%./nt;
%m.lp_lrms=lp_lrms./nt;
m.hp_0=hp_0;%./nt;
%m.hp_lrms=hp_lrms./nt;


m.volfrac_lp_0   =volfrac_lp_0;%./nt;
%m.volfrac_lp_lrms=volfrac_lp_lrms./nt;
m.volfrac_hp_0   =volfrac_hp_0;%./nt;
%m.volfrac_hp_lrms=volfrac_hp_lrms./nt;
