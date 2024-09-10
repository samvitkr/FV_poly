jcond=188;

fil=sprintf('../data/lsefilter_j_%03d.mat',jcond);
mfil=matfile(fil);

fn=sprintf('../data/velgradfield_lseQ4_j_%03d.mat',jcond);
m=matfile(fn);

fnd=sprintf('../data/velgradfield_lseQ4_dfil_j_%03d.mat',jcond);
md=matfile(fnd,'Writable',true);

md.u=ifft2( fft2(m.u).*mfil.dfil ,'symmetric' );
md.v=ifft2( fft2(m.v).*mfil.dfil ,'symmetric' );
md.w=ifft2( fft2(m.w).*mfil.dfil ,'symmetric' );

md.dudx=ifft2( fft2(m.dudx).*mfil.dfil ,'symmetric' );
md.dvdx=ifft2( fft2(m.dvdx).*mfil.dfil ,'symmetric' );
md.dwdx=ifft2( fft2(m.dwdx).*mfil.dfil ,'symmetric' );

md.dudy=ifft2( fft2(m.dudy).*mfil.dfil ,'symmetric' );
md.dvdy=ifft2( fft2(m.dvdy).*mfil.dfil ,'symmetric' );
md.dwdy=ifft2( fft2(m.dwdy).*mfil.dfil ,'symmetric' );

md.dudz=ifft2( fft2(m.dudz).*mfil.dfil ,'symmetric' );
md.dvdz=ifft2( fft2(m.dvdz).*mfil.dfil ,'symmetric' );
md.dwdz=ifft2( fft2(m.dwdz).*mfil.dfil ,'symmetric' );

fnu=sprintf('../data/velgradfield_lseQ4_ufil_j_%03d.mat',jcond);
mu=matfile(fnu,'Writable',true);
mu.u=ifft2( fft2(m.u).*mfil.ufil ,'symmetric' );
mu.v=ifft2( fft2(m.v).*mfil.ufil ,'symmetric' );
mu.w=ifft2( fft2(m.w).*mfil.ufil ,'symmetric' );

mu.dudx=ifft2( fft2(m.dudx).*mfil.ufil ,'symmetric' );
mu.dvdx=ifft2( fft2(m.dvdx).*mfil.ufil ,'symmetric' );
mu.dwdx=ifft2( fft2(m.dwdx).*mfil.ufil ,'symmetric' );

mu.dudy=ifft2( fft2(m.dudy).*mfil.ufil ,'symmetric' );
mu.dvdy=ifft2( fft2(m.dvdy).*mfil.ufil ,'symmetric' );
mu.dwdy=ifft2( fft2(m.dwdy).*mfil.ufil ,'symmetric' );                                                                               

mu.dudz=ifft2( fft2(m.dudz).*mfil.ufil ,'symmetric' );
mu.dvdz=ifft2( fft2(m.dvdz).*mfil.ufil ,'symmetric' );
mu.dwdz=ifft2( fft2(m.dwdz).*mfil.ufil ,'symmetric' );


