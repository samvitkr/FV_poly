%time=52000;

load('dpdx.mat')
%load('volc.mat');
ut=mean(ut_ts);
%be=[-inf,-100:1:100,inf]/500;
%be= (ut^2)*[-inf,-100:0.5:100,inf];
%be= (ut^2)*[-inf,-300:0.5:300,inf];
%be= (ut^2)*[-inf,-2000:5:2000,inf];
be= (ut^2)*[-inf,-5000:10:5000,inf];
nbins=length(be)-1;

m=matfile('hist_l_mean.mat','Writable',true)
vol=zeros(1,nbins);
contrip=zeros(1,nbins);
count=0;
for time=0:1000:40000;
	count=count+1;
	fh=sprintf("hist_l_%07d.mat",time)
	mh=matfile(fh)
	%mh.BinEdgeslambda=be;
	vol=vol+mh.vol;
	contrip=contrip+mh.contrip;
end
m.vol=vol./count;
m.contrip=contrip./count;
m.BinEdgeslambda=be;
m.ut=ut;
