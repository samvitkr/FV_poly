%time=52000;

load('dpdx.mat')
load('volc.mat');
%ut=ut_ts(time+1);
ut=mean(ut_ts);
be=[-inf,-100:1:100,inf];
nbins=length(be)-1;

for time=40000:1000:80000
	ft=sprintf("transferfields_%07d.mat",time)
	m=matfile(ft)
	%fp=sprintf("profiles_%07d.mat",time)
	%mp=matfile(fp,'Writable',true)
	%load("volc.mat");
	%nbins=1000;
	%mp.vozm=squeeze(mean(mean(m.voz,1),2));
	%mp.woym=squeeze(mean(mean(m.woy,1),2));
	%mp.viscm=squeeze(mean(mean(m.visc,1),2));
	%mp.polym=squeeze(mean(mean(m.poly,1),2));
	syz=m.voz-m.woy+m.visc+m.poly;
	syzc=0.5*(syz(:,:,1:end-1)+syz(:,:,2:end));
	syzcp=syzc./(-ut^2);
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
		i=logical((syzcp>e1).*(syzcp<e2));
		volh(n)=sum(volc(i));
		contri(n)=sum(integ(i));
	end
	fh=sprintf("hist_%07d.mat",time)
	mh=matfile(fh,'Writable',true)
	mh.BinEdges=be;
	mh.vol=volh;
	mh.contrip=contri;
end
