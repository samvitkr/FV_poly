%time=52000;

load('dpdx.mat')
%load('volc.mat');
ut=mean(ut_ts);
%be=[-inf,-100:1:100,inf];
be=(ut^2)*[-inf,-100:0.5:100,inf];
nbins=length(be)-1;

m=matfile('hist_wi_6p67_l_mean.mat','Writable',true)
vol=zeros(1,nbins);
contrip=zeros(1,nbins);
count=0;
for time=0:1000:60000;
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
