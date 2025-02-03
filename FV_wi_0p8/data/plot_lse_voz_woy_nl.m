close all
clear
jcond=156;
load('../data/ygrid.mat')
mp=matfile("../data/mean_profiles.mat")
ft=sprintf('../data/velfield_lse_voz_uufilter_j_%03d.mat',jcond);
m=matfile(ft)
f1=sprintf("lse_voz_uufilter_profiles_j_%03d.fig",jcond)
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

woy2m=squeeze(mean(woy2,[1 2]));
woy4m=squeeze(mean(woy4,[1 2]));
voz2m=squeeze(mean(voz2,[1 2]));
voz4m=squeeze(mean(voz4,[1 2]));
nl2m=squeeze(mean(nl2,[1 2]));
nl4m=squeeze(mean(nl4,[1 2]));

voz2tm=squeeze(mean(voz2t,[1 2]));
voz4tm=squeeze(mean(voz4t,[1 2]));
nl2tm=squeeze(mean(nl2t,[1 2]));
nl4tm=squeeze(mean(nl4t,[1 2]));

y=yCheb+1;
h1=figure
subplot(2,3,1)
hold on
plot(-voz2m,y,'--b','LineWidth',1.5)
plot(-voz4m,y,'--r','LineWidth',1.5)
hold off
ylim([0 1])
yline(y(156))
xlim([-1e-5 1e-5])
xline(0)
xlabel("-v\omega_z' ")
subplot(2,3,2)
hold on
plot(woy2m,y,':b','LineWidth',1.5)
plot(woy4m,y,':r','LineWidth',1.5)
hold off
ylim([0 1])
yline(y(156))
xlim([-1e-5 1e-5])
xline(0)
xlabel('w\omega_y')

subplot(2,3,3)
hold on
plot(-nl2m,y,'-b','LineWidth',1.5)
plot(-nl4m,y,'-r','LineWidth',1.5)
hold off
ylim([0 1])
yline(y(156))
xlim([-1e-5 1e-5])
xline(0)
xlabel("-(v\omega_z'-w\omega_y)")

subplot(2,3,4)
hold on
plot(-voz2tm,y,'-.b','LineWidth',1.5)
plot(-voz4tm,y,'-.r','LineWidth',1.5)
hold off
ylim([0 1])
yline(y(156))
xlim([-1e-5 1e-5])
xline(0)
xlabel("-v\omega_z ")

subplot(2,3,5)
hold on
plot(-nl2tm,y,'-b','LineWidth',1.5)
plot(-nl4tm,y,'-r','LineWidth',1.5)
hold off
ylim([0 1])
yline(y(156))
xlim([-1e-5 1e-5])
xline(0)
xlabel("-(v\omega_z-w\omega_y)")
saveas(h1,f1)