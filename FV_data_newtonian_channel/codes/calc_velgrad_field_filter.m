tstart=40000;
tend=40000;
tstep=1000;
Nx=512;
Ny=220;
Nz=384;
mf=matfile('../data/filter.mat')

for time=tstart:tstep:tend
        %fvg=sprintf('velgrad_transfer_flp_%07d.mat',time);
%       fvg=sprintf('velgrad_transfer_fhp_%07d.mat',time);
	fv=sprintf('../data/velfields_%07d.mat',time);
        fvg=sprintf('../data/velgrad_%07d.mat',time);
	ft=sprintf('../data/transferfields_%07d.mat',time);
	fgd=sprintf('../data/velgrad_transfer_dfil_%07d.mat',time);
	fgu=sprintf('../data/velgrad_transfer_ufil_%07d.mat',time);
        %        fo=sprintf('vortfields_%07d.mat',time);
	mv=matfile(fv)
	mvg=matfile(fvg)
	mt=matfile(ft)
	mgd=matfile(fgd,'Writable',true);
	mgu=matfile(fgu,'Writable',true);

	mgd.u=ifft2(fft2(mv.ufield).*mf.fulldfil,'symmetric');
	mgd.v=ifft2(fft2(mv.vfield).*mf.fulldfil,'symmetric');
	mgd.w=ifft2(fft2(mv.wfield).*mf.fulldfil,'symmetric');
	mgd.dudx=ifft2(fft2(mvg.dudx).*mf.fulldfil,'symmetric');
	mgd.dvdx=ifft2(fft2(mvg.dvdx).*mf.fulldfil,'symmetric');
	mgd.dwdx=ifft2(fft2(mvg.dwdx).*mf.fulldfil,'symmetric');
	mgd.dudy=ifft2(fft2(mvg.dudy).*mf.fulldfil,'symmetric');
	mgd.dvdy=ifft2(fft2(mvg.dvdy).*mf.fulldfil,'symmetric');
	mgd.dwdy=ifft2(fft2(mvg.dwdy).*mf.fulldfil,'symmetric');
	mgd.dudz=ifft2(fft2(mvg.dudz).*mf.fulldfil,'symmetric');
	mgd.dvdz=ifft2(fft2(mvg.dvdz).*mf.fulldfil,'symmetric');
	mgd.dwdz=ifft2(fft2(mvg.dwdz).*mf.fulldfil,'symmetric');
	mgd.visc=ifft2(fft2(mt.visc).*mf.fulldfil,'symmetric');


	mgu.u=ifft2(fft2(mv.ufield).*mf.fullufil,'symmetric');
        mgu.v=ifft2(fft2(mv.vfield).*mf.fullufil,'symmetric');
        mgu.w=ifft2(fft2(mv.wfield).*mf.fullufil,'symmetric');
        mgu.dudx=ifft2(fft2(mvg.dudx).*mf.fullufil,'symmetric');
        mgu.dvdx=ifft2(fft2(mvg.dvdx).*mf.fullufil,'symmetric');
        mgu.dwdx=ifft2(fft2(mvg.dwdx).*mf.fullufil,'symmetric');
        mgu.dudy=ifft2(fft2(mvg.dudy).*mf.fullufil,'symmetric');
        mgu.dvdy=ifft2(fft2(mvg.dvdy).*mf.fullufil,'symmetric');
        mgu.dwdy=ifft2(fft2(mvg.dwdy).*mf.fullufil,'symmetric');
        mgu.dudz=ifft2(fft2(mvg.dudz).*mf.fullufil,'symmetric');
        mgu.dvdz=ifft2(fft2(mvg.dvdz).*mf.fullufil,'symmetric');
        mgu.dwdz=ifft2(fft2(mvg.dwdz).*mf.fullufil,'symmetric');
        mgu.visc=ifft2(fft2(mt.visc).*mf.fullufil,'symmetric');
	
end
