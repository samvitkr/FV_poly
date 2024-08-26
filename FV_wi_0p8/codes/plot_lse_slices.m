close all
clear
jcond=156;
load('../data/ygrid.mat')
mp=matfile("../data/mean_profiles.mat")
%ft=sprintf('../data/velfield_lse_voz_ddfilter_j_%03d.mat',jcond)
ft=sprintf('../data/velfield_lse_j_%03d.mat',jcond)
m=matfile(ft)
f1=sprintf("lse_voz_ddfilter_slice_j_%03d.fig",jcond)
ozm=reshape(-mp.dUdy,[1 1 220]);
woy2=(m.wQ2).*(m.Q2omegay);
woy4=(m.wQ4).*(m.Q4omegay);
voz2=(m.vQ2).*(m.Q2omegaz);
voz4=(m.vQ4).*(m.Q4omegaz);
voz2t=(m.vQ2).*(m.Q2omegaz+ozm);
voz4t=(m.vQ4).*(m.Q4omegaz+ozm);

nl2=voz2-woy2;
nl4=voz4-woy4;
nl2t=voz2t-woy2;
nl4t=voz4t-woy4;
%%
close all

%%
h1=figure
subplot(2,2,1)
pcolor(m.X(:,:,jcond),m.Z(:,:,jcond),-voz2t(:,:,jcond))
shading flat
h=colorbar
ylabel(h,'-v\omega_z')
clim([-1e-3 1e-3])

subplot(2,2,2)
pcolor(m.X(:,:,jcond),m.Z(:,:,jcond),-voz4t(:,:,jcond))
shading flat
clim([-1e-3 1e-3])
h=colorbar
ylabel(h,'-v\omega_z')
colormap jet

subplot(2,2,3)
pcolor(m.X(:,:,jcond),m.Z(:,:,jcond),woy2(:,:,jcond))
shading flat
h=colorbar
ylabel(h,'w\omega_y')
clim([-5e-4 5e-4])

subplot(2,2,4)
pcolor(m.X(:,:,jcond),m.Z(:,:,jcond),woy4(:,:,jcond))
shading flat
clim([-5e-4 5e-4])
h=colorbar
ylabel(h,'w\omega_y')
colormap jet
%saveas(h1,f1)
%%
close all
v=m.vQ2(:,:,jcond);
oz=m.Q2omegaz(:,:,jcond);
w=m.wQ2(:,:,jcond);
oy=m.Q2omegay(:,:,jcond);
subplot(1,2,1)
phiv=fft2(v);
phioz=fft2(oz);
phivoz=phiv.*conj(phioz);
phioy=fft2(oy);
phiw=fft2(w);
phiwoy=phiw.*conj(phioy);
pcolor(real(phivoz(1:192,1:256)-phiwoy(1:192,1:256)))
shading flat
clim([-1000 1000])
colormap jet
v=m.vQ4(:,:,jcond);
oz=m.Q4omegaz(:,:,jcond);
w=m.wQ4(:,:,jcond);
oy=m.Q4omegay(:,:,jcond);
subplot(1,2,2)
phiv=fft2(v);
phioz=fft2(oz);
phivoz=phiv.*conj(phioz);
phioy=fft2(oy);
phiw=fft2(w);
phiwoy=phiw.*conj(phioy);
pcolor(real(phivoz(1:192,1:256)-phiwoy(1:192,1:256)))
shading flat
clim([-1000 1000])
colormap jet