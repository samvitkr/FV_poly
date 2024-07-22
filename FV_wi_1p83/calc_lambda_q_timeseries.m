tstart=40000;
tend=130000;
tstep=1000;
Nx=512;
Ny=220;
Nz=384;
Lx=4*pi;
Lz=2*pi;
load('ygrid.mat')
load('dpdx.mat')
nt=(tend-tstart)/tstep+1;
lambda=zeros(nt,1);
Q=zeros(nt,1);
Qp=zeros(nt,1);
SYZ=zeros(nt,1);
NL=zeros(nt,1);
NLABS=zeros(nt,1);
NLABSF=zeros(nt,1);
POLY=zeros(nt,1);
VISC=zeros(nt,1);
POLYABS=zeros(nt,1);
VISCABS=zeros(nt,1);
ts=0;
for time=tstart:tstep:tend
	
        f=sprintf("lambda_%07d.mat",time)
        ml=matfile(f);
	
	ft=sprintf("transferfields_%07d.mat",time)
	mt=matfile(ft);
	nl=mt.voz-mt.woy;
	poly=mt.poly;
	visc=mt.visc;
	syz=nl+poly+visc;

	nlm=squeeze(mean(nl,[1,2]));
        polym=squeeze(mean(poly,[1,2]));
        viscm=squeeze(mean(visc,[1,2]));
	syzm=squeeze(mean(syz,[1,2]));
	nlabs=abs(nlm);
	nlabsf=squeeze(mean(abs(nl),[1,2]));
	polyabsm=squeeze(mean(  abs(poly),[1,2]));
        viscabsm=squeeze(mean(  abs(visc),[1,2]));


	l=ml.lambda2;
	q=ml.Q;
	qpos=ml.Q>0;
        qp=qpos.*(ml.Q);
        ln=ml.lambda2;
	lm=squeeze(mean(l,[1,2]));
	qm=squeeze(mean(q,[1,2]));
	qpm=squeeze(mean(qp,[1,2]));
	y=flip(yCheb);
	li=0.5*trapz(y,lm);
	qi=0.5*trapz(y,qm);
	qpi=0.5*trapz(y,qpm);
	nli=0.5*trapz(y,nlm);
	nlabsi=0.5*trapz(y,nlabs);
        nlabsif=0.5*trapz(y,nlabsf);	
        polyi=0.5*trapz(y,polym);
        visci=0.5*trapz(y,viscm);
	polyabsi=0.5*trapz(y,polyabsm);
        viscabsi=0.5*trapz(y,viscabsm);
	syzi=0.5*trapz(y,syzm);
	ts=ts+1;
	lambda(ts)=li;
	Q(ts)=qi;
	Qp(ts)=qpi;
	NL(ts)=nli;
	NLABS(ts)=nlabsi;
	NLABSF(ts)=nlabsif;
	SYZ(ts)=syzi;
	POLY(ts)=polyi;
	VISC(ts)=visci;
	POLYABS(ts)=polyabsi;
        VISCABS(ts)=viscabsi;
end
m=matfile('lambda_Q_timeseries.mat','Writable',true)
m.lambda=lambda;
m.Q=Q;
m.Qp=Qp;
m.nl=NL;
m.nlabs=NLABS;
m.nlabsf=NLABSF;
m.syz=SYZ;
m.poly=POLY;
m.visc=VISC;
m.polyabs=POLYABS;
m.viscabs=VISCABS;
