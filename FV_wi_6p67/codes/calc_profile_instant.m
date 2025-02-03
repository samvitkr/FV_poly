for time=000:1000:120000
	ft=sprintf("transferfields_%07d.mat",time)
	m=matfile(ft)
	fv=sprintf("velfields_%07d.mat",time)
	mv=matfile(fv);
	fp=sprintf("./profiles/profiles_%07d.mat",time)
	mp=matfile(fp,'Writable',true)
	mp.vozm=squeeze(mean(mean(m.voz,1),2));
	mp.woym=squeeze(mean(mean(m.woy,1),2));
	mp.viscm=squeeze(mean(mean(m.visc,1),2));
	mp.polym=squeeze(mean(mean(m.poly,1),2));
	mp.Um=squeeze(mean(mean(mv.ufield,1),2));
%	mp.lm=exp(squeeze(mean(mean(  log(m.l) ,1),2)));
end
