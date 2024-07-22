re=4667;
ret=278;
ut=ret/re;
dnu=1/ret;
load('ygrid.mat')
m=matfile('transferfields_mean.mat')

vozm=squeeze(mean(mean(m.voz,1),2));
woym=squeeze(mean(mean(m.woy,1),2));
viscm=squeeze(mean(mean(m.visc,1),2));
%polym=squeeze(mean(mean(m.poly,1),2));
Um=squeeze(mean(mean(m.U,1),2));
uu=squeeze(mean(mean(m.uu,1),2))-Um.^2;
vv=squeeze(mean(mean(m.vv,1),2));
ww=squeeze(mean(mean(m.ww,1),2));
uv=squeeze(mean(mean(m.uv,1),2));
uw=squeeze(mean(mean(m.uw,1),2));
vw=squeeze(mean(mean(m.vw,1),2));

m2=matfile('mean_profiles.mat','Writable',true)
m2.vozm=vozm;
m2.woym=woym;
 %m2.polym=polym;
m2.viscm=viscm;
m2.Um=Um;
m2.yCheb=yCheb;
%m2.ret=ret;
%m2.re=re;
m2.uu=uu;
m2.vv=vv;
m2.ww=ww;
m2.uv=uv;
m2.uw=uw;
m2.vw=vw;


%x1=180;
%y1=180;
%x2=250;
%y2=300;
%
%h1=figure('OuterPosition',...
%    [x1 y1 x2 y2]);
%hold on
%plot(vozm./(-ut^2),yCheb,'--b','LineWidth',1.5)
%plot(-woym./(-ut^2),yCheb,':b','LineWidth',1.5)
%plot((vozm-woym)./(-ut^2),yCheb,'-b','LineWidth',1.5)
%hold off
%ylabel('y/H')
%xlabel('$\times (-u_\tau^2)/H$','interpreter','latex')
%xline(1)
%pbaspect([1 1.5 1])
%
%legend('v\omega_z','-w\omega_y','v\omega_z-w\omega_y','location','west')
%legend boxoff
%set(gca,'FontSize',10)
%grid on
%print(h1,'mean_conv','-dpng');
%
%close all
%
%
%h2=figure('OuterPosition',...
%    [x1 y1 x2 y2]);
%hold on
%plot(viscm./(-ut^2),yCheb,'-r','LineWidth',1.5)
%%plot(polym./(-ut^2),yCheb,'-m','LineWidth',1.5)
%hold off
%ylabel('y/H')
%xlabel('$\times (-u_\tau^2)/H$','interpreter','latex')
%xline(1)
%pbaspect([1 1.5 1])
%
%legend('visc','poly','location','east')
%legend boxoff
%set(gca,'FontSize',10)
%grid on
%print(h2,'mean_visc_poly','-dpng');
%
%close all
%
%h3=figure('OuterPosition',...
%    [x1 y1 x2 y2]);
%hold on
%plot( (vozm-woym+viscm)./(-ut^2),yCheb,'-k','LineWidth',1.5)
%hold off
%
%xlim([0.5 1.5])
%set(gca,'XMinorTick','on')
%ylabel('y/H')
%xlabel('$\times (-u_\tau^2)/H$','interpreter','latex')
%xline(1)
%pbaspect([1 1.5 1])
%
%legend('\Sigma_{yz}','location','east')
%legend boxoff
%set(gca,'FontSize',10)
%grid on
%print(h3,'mean_tot','-dpng');











