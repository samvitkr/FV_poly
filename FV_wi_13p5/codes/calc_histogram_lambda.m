%time=52000;

load('../data/dpdx.mat')
load('../data/volc.mat');
%ut=ut_ts(time+1);
ut=mean(ut_ts);
%be=[-inf,-100:1:100,inf]/500;
%be= (ut^2)*[-inf,-100:0.5:100,inf];
be=[-inf,-0.2:0.01:0.2,inf];
nbins=length(be)-1;
ts1=1000*[92,158,234,256,316]/2;
ts2=1000*[104,178,244,264,330]/2;
ts=[ts1,ts2];
ts=sort(ts);
nts=length(ts);
for i=1:nts;
	time=ts(i);
	ft=sprintf("../data/transferfields_%07d.mat",time);
	m=matfile(ft);
	fl=sprintf("../data/lambda_%07d.mat",time);
	ml=matfile(fl);
	l=ml.lambda2;
	% %load("volc.mat");
	% %nbins=1000;
	%mp.vozm=squeeze(mean(mean(m.voz,1),2));
	%mp.woym=squeeze(mean(mean(m.woy,1),2));
	%mp.viscm=squeeze(mean(mean(m.visc,1),2));
	%mp.polym=squeeze(mean(mean(m.poly,1),2));
	syz=m.voz-m.woy+m.visc+m.poly;
	syzc=0.5*(syz(:,:,1:end-1)+syz(:,:,2:end));
	syzcp=syzc./(-ut^2);
	lc=0.5*(l(:,:,1:end-1)+l(:,:,2:end));
	
	%size(l)
	%size(lc)
	%h=histogram(syzc,nbins);
	%be=h.BinEdges;
	%volh=h.Values*0;
	%be=[-inf,-100:1:100,inf];
	%nbins=length(be)-1;
	integ=syzcp.*volc;
	for n=1:nbins
		n
		e1=be(n);
		e2=be(n+1);
		i=logical((lc>e1).*(lc<e2));
		volh(n)=sum(volc(i));
		contri(n)=sum(integ(i));
	end
	fh=sprintf("../data/hist_l_%07d.mat",time)
	mh=matfile(fh,'Writable',true)
	mh.BinEdgeslambda=be;
	mh.vol=volh;
	mh.contrip=contri;
end
