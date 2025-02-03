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
        %l=ml.Q;
	flag=ml.Q>0.1;
	voz=squeeze(mean(m.voz.*flag,[1 2]));
	woy=squeeze(mean(m.woy.*flag,[1 2]));
	visc=squeeze(mean(m.visc.*flag,[1 2]));
	poly=squeeze(mean(m.poly.*flag,[1 2]));
	volfrac=squeeze(mean(flag,[1 2]));
	fh=sprintf("../data/profile_q0p1_%07d.mat",time)
        mh=matfile(fh,'Writable',true)
	mh.voz=voz;
	mh.woy=woy;
	mh.visc=visc;
	mh.poly=poly;
	mh.volfrac=volfrac;

end
