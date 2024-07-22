tstart=10000;
tend=50000;
tstep=1000;
Nx=512;
Ny=220;
Nz=384;
nt=(tend-tstart)/tstep+1;
visc_0_pi_top=zeros(nt,1);
visc_pi_2pi_top=zeros(nt,1);
visc_0_pi_bot=zeros(nt,1);
visc_pi_2pi_bot=zeros(nt,1);


poly_0_pi_top=zeros(nt,1);
poly_pi_2pi_top=zeros(nt,1);
poly_0_pi_bot=zeros(nt,1);
poly_pi_2pi_bot=zeros(nt,1);



mn=matfile('transferfields_quad_l_mean_opt.mat','Writable',true)
ks=mn.ks;

ts=0;
for time=tstart:tstep:tend
        ts=ts+1;
	ft=sprintf("sf_%07d.mat",time)
        mt=matfile(ft);
	km=ks(ts);
	%km=0;
	visc_bot=circshift(mt.visc_bot,km,1);
        visc_top=circshift(mt.visc_top,km,1);
	poly_bot=circshift(mt.poly_bot,km,1);
	poly_top=circshift(mt.poly_top,km,1);

	visc_0_pi_top(ts)=0.25*mean( visc_top(1:Nz/2,:),'all' );
	visc_pi_2pi_top(ts)=0.25*mean( visc_top(1+Nz/2:end,:),'all' ) ;
	visc_0_pi_bot(ts)=0.25*mean( visc_bot(1:Nz/2,:),'all' );
	visc_pi_2pi_bot(ts)=0.25*mean( visc_bot(1+Nz/2:end,:),'all' );


	poly_0_pi_top(ts)=0.25*mean( poly_top(1:Nz/2,:),'all' );
        poly_pi_2pi_top(ts)=0.25*mean( poly_top(1+Nz/2:end,:),'all' ) ;
        poly_0_pi_bot(ts)=0.25*mean( poly_bot(1:Nz/2,:),'all' );
        poly_pi_2pi_bot(ts)=0.25*mean( poly_bot(1+Nz/2:end,:),'all' );
end

mn.sf_0_pi_top=  visc_0_pi_top  +poly_0_pi_top;
mn.sf_pi_2pi_top=visc_pi_2pi_top+poly_pi_2pi_top;
mn.sf_0_pi_bot=  visc_0_pi_bot  +poly_0_pi_bot;
mn.sf_pi_2pi_bot=visc_pi_2pi_bot+poly_pi_2pi_bot ;


mn.visc_0_pi_top=  visc_0_pi_top  ;
mn.visc_pi_2pi_top=visc_pi_2pi_top;
mn.visc_0_pi_bot=  visc_0_pi_bot  ;
mn.visc_pi_2pi_bot=visc_pi_2pi_bot;

mn.poly_0_pi_top=  poly_0_pi_top;
mn.poly_pi_2pi_top=poly_pi_2pi_top;
mn.poly_0_pi_bot=  poly_0_pi_bot;
mn.poly_pi_2pi_bot=poly_pi_2pi_bot ;
