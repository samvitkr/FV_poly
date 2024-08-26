close all
clear
jcond=156;
lt=5e-3;
%fth=0.1;
%ft=sprintf("velfield_lse_voz_j_%d.mat",jcond)
%ft=sprintf("velfield_lse_vwoy_j_%d.mat",jcond)
load('../data/ygrid.mat')
ft=sprintf('../data/velfield_lse_ddfilter_j_%03d.mat',jcond);
m=matfile(ft)
ftu=sprintf('../data/velfield_lse_uufilter_j_%03d.mat',jcond);
mu=matfile(ftu)


q4=m.Q4Q;
qu4=mu.Q4Q;
q2=m.Q2Q;
qu2=mu.Q2Q;

voz2=(m.vQ2).*(m.Q2omegaz);
woy2=(m.wQ2).*(m.Q2omegay);
voz4=(m.vQ4).*(m.Q4omegaz);
woy4=(m.wQ4).*(m.Q4omegay);
nl2=voz2-woy2;
nl4=voz4-woy4;

voz2u=(mu.vQ2).*(mu.Q2omegaz);
woy2u=(mu.wQ2).*(mu.Q2omegay);
voz4u=(mu.vQ4).*(mu.Q4omegaz);
woy4u=(mu.wQ4).*(mu.Q4omegay);
nl2u=voz2u-woy2u;
nl4u=voz4u-woy4u;

nl=nl4;
nlu=nl4u;
q=q4;
qu=qu4;

x1=150;
y1=150;
x2=1000;
y2=350;
h1=figure('OuterPosition',...
    [x1 y1 x2 y2]);
subplot(1,2,1)
%isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
%permute(q,[2 1 3]), lt, permute((1e+3)*m.fxQ2,[2 1 3]))

isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
permute(q,[2 1 3]), lt, permute(-nl,[2 1 3]))

%isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
%permute((1e+3)*m.fxQ2,[2 1 3]), -fth)
%isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
%permute((1e+3)*m.fxQ2,[2 1 3]), fth)
%clim([-1 1])
colormap jet
colorbar
lightangle(-45,-90)
axis equal
ylim([-1 1])
xlim([-0.5 0.5])
% zlim([0 0.4])
view(45,45)
xlabel('z')
ylabel('x')
zlabel('y')
grid on

subplot(1,2,2)
%isosurface( permute(mu.Z,[2 1 3]), permute(mu.X,[2 1 3]), permute(mu.Y,[2 1 3]),...
%permute(qu,[2 1 3]), lt, permute((1e+3)*mu.fxQ2,[2 1 3])) 

isosurface( permute(mu.Z,[2 1 3]), permute(mu.X,[2 1 3]), permute(mu.Y,[2 1 3]),...
permute(qu,[2 1 3]), lt, permute(-nlu,[2 1 3])) 

%isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
%permute((1e+3)*m.fyQ2,[2 1 3]), -fth)
%isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
%permute((1e+3)*m.fyQ2,[2 1 3]), fth)
%clim([-1 1])
colormap jet
colorbar
lightangle(-45,-90)
axis equal
ylim([-1 1])
xlim([-0.5 0.5])
% zlim([0 0.4])
view(45,45)
xlabel('z')
ylabel('x')
zlabel('y')
grid on


f1=sprintf("iso_lse_vel_nl_Q5e-3_filter_j2_%03d.fig",jcond)
saveas(h1,f1)
fl=sprintf("lse_profiles_j_%03d.mat",jcond)
mfl=matfile(fl,'Writable',true)
mfl.lt=lt;

mfl.voz2=squeeze(mean(voz2,[1 2]));
mfl.voz4=squeeze(mean(voz4,[1 2]));
mfl.woy2=squeeze(mean(woy2,[1 2]));
mfl.woy4=squeeze(mean(woy4,[1 2]));

mfl.voz2u=squeeze(mean(voz2u,[1 2]));
mfl.voz4u=squeeze(mean(voz4u,[1 2]));
mfl.woy2u=squeeze(mean(woy2u,[1 2]));
mfl.woy4u=squeeze(mean(woy4u,[1 2]));

mfl.voz2q=squeeze(mean(voz2.*(q2>lt),[1 2]));
mfl.voz4q=squeeze(mean(voz4.*(q4>lt),[1 2]));
mfl.woy2q=squeeze(mean(woy2.*(q2>lt),[1 2]));
mfl.woy4q=squeeze(mean(woy4.*(q4>lt),[1 2]));

mfl.voz2uq=squeeze(mean(voz2u.*(qu2>lt),[1 2]));
mfl.voz4uq=squeeze(mean(voz4u.*(qu4>lt),[1 2]));
mfl.woy2uq=squeeze(mean(woy2u.*(qu2>lt),[1 2]));
mfl.woy4uq=squeeze(mean(woy4u.*(qu4>lt),[1 2]));

mfl.fx2=squeeze(mean(m.fxQ2,[1 2]));
mfl.fx4=squeeze(mean(m.fxQ4,[1 2]));
mfl.fx2u=squeeze(mean(mu.fxQ2,[1 2]));
mfl.fx4u=squeeze(mean(mu.fxQ4,[1 2]));

mfl.fx2q=squeeze(mean(m.fxQ2.*(q2>lt),[1 2]));
mfl.fx4q=squeeze(mean(m.fxQ4.*(q4>lt),[1 2]));
mfl.fx2uq=squeeze(mean(mu.fxQ2.*(qu2>lt),[1 2]));
mfl.fx4uq=squeeze(mean(mu.fxQ4.*(qu4>lt),[1 2]));

% close all
%%
% q=m.Q4Q;
% qu=mu.Q4Q;
% %mq1=mean(abs(q),'all')
% %mq2=mean(abs(qu),'all')
% nl4=(m.vQ4).*(m.Q4omegaz)-(m.wQ4).*(m.Q4omegay);
% nlu4=(mu.vQ4).*(mu.Q4omegaz)-(mu.wQ4).*(mu.Q4omegay);
% nlq4=(nl4).*(q>lt);
% nlqu4=(nlu4).*(qu>lt);
% nlm4=squeeze(mean(nl4,[1 2]));
% nlum4=squeeze(mean(nlu4,[1 2]));
% nlqm4=squeeze(mean(nl4,[1 2]));
% nlqum4=squeeze(mean(nlu4,[1 2]));
% 
% q=m.Q2Q;
% qu=mu.Q2Q;
% nl2=(m.vQ2).*(m.Q2omegaz)-(m.wQ2).*(m.Q2omegay);
% nlu2=(mu.vQ2).*(mu.Q2omegaz)-(mu.wQ2).*(mu.Q2omegay);
% nlq2=(nl2).*(q>lt);
% nlqu2=(nlu2).*(qu>lt);
% nlm2=squeeze(mean(nl2,[1 2]));
% nlum2=squeeze(mean(nlu2,[1 2]));
% nlqm2=squeeze(mean(nl2,[1 2]));
% nlqum2=squeeze(mean(nlu2,[1 2]));


%%
% figure
% hold on
% plot(-nlqm2,yCheb+1,'-b','LineWidth',1.5)
% plot(-nlqum2,yCheb+1,'-r','LineWidth',1.5)
% plot(-nlqm4,yCheb+1,'--b','LineWidth',1.5)
% plot(-nlqum4,yCheb+1,'--r','LineWidth',1.5)
% hold off
% yline(yCheb(jcond)+1)
% ylim([0 1])
% grid on
%%
% h2=figure('OuterPosition',...
%     [x1 y1 x2 y2]);
% subplot(1,3,1)
% isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
% permute(m.Q4lambda2,[2 1 3]), lt, permute(260*m.uQ4, [2 1 3]))
% % isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
% % permute((1e+3)*m.fxQ2,[2 1 3]), -fth)
% % isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
% % permute((1e+3)*m.fxQ2,[2 1 3]), fth)
% clim([-1 1])
% colormap jet
% colorbar
% lightangle(-45,-90)
% axis equal
% ylim([-1 1])
% xlim([-0.5 0.5])
% % zlim([0 0.4])
% view(45,45)
% xlabel('z')
% ylabel('x')
% zlabel('y')
% grid on
% 
% subplot(1,3,2)
% isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
% permute(m.Q4lambda2,[2 1 3]), lt, permute(260*m.vQ4, [2 1 3]))
% % isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
% % permute((1e+3)*m.fxQ2,[2 1 3]), -fth)
% % isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
% % permute((1e+3)*m.fxQ2,[2 1 3]), fth)
% clim([-1 1])
% colormap jet
% colorbar
% lightangle(-45,-90)
% axis equal
% ylim([-1 1])
% xlim([-0.5 0.5])
% % zlim([0 0.4])
% view(45,45)
% xlabel('z')
% ylabel('x')
% zlabel('y')
% grid on
% 
% subplot(1,3,3)
% isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
% permute(m.Q4lambda2,[2 1 3]), lt, permute(260*m.wQ4, [2 1 3]))
% % isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
% % permute((1e+3)*m.fxQ2,[2 1 3]), -fth)
% % isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
% % permute((1e+3)*m.fxQ2,[2 1 3]), fth)
% clim([-1 1])
% colormap jet
% colorbar
% lightangle(-45,-90)
% axis equal
% ylim([-1 1])
% xlim([-0.5 0.5])
% % zlim([0 0.4])
% view(45,45)
% xlabel('z')
% ylabel('x')
% zlabel('y')
% grid on
% f1=sprintf("iso_lse_vwoy2_j_%03d.fig",jcond)
% saveas(h2,f1)
%%
% subplot(1,2,2)
% isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
% permute(m.Q4lambda2,[2 1 3]), lt, permute(260*m.uQ4, [2 1 3]))
% isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
% permute((1e+3)*m.fxQ4,[2 1 3]), -fth)
% isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
% permute((1e+3)*m.fxQ4,[2 1 3]), fth)
% clim([-1.5 1.5])
% colormap jet
% colorbar
% lightangle(-45,-90)
% axis equal
% % xlim([-0.5 0.5])
% % ylim([-0.5 0.5])
% % zlim([0 0.4])
% view(45,45)
% xlabel('z')
% ylabel('x')
% zlabel('y')
% %  f1=sprintf("iso_lse_voz_j_%03d.fig",jcond)
% %  saveas(h1,f1)
 %%
% 
%  ft=sprintf("velfield_lse_woy_j_%d.mat",jcond)
%  m=matfile(ft)
% 
%  h1=figure('OuterPosition',...
%     [x1 y1 x2 y2]);
% subplot(1,2,1)
% isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
% permute(m.Q2lambda2,[2 1 3]), lt, permute(260*m.wQ2, [2 1 3]))
% isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
% permute((1e+3)*m.fxQ2,[2 1 3]), -fth)
% isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
% permute((1e+3)*m.fxQ2,[2 1 3]), fth)
% clim([-1.5 1.5])
% colormap jet
% colorbar
% lightangle(-45,-90)
% axis equal
% xlim([-0.5 0.5])
% ylim([-0.5 0.5])
% zlim([0 0.4])
% view(45,45)
% xlabel('z')
% ylabel('x')
% zlabel('y')
% 
% subplot(1,2,2)
% isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
% permute(m.Q4lambda2,[2 1 3]), lt, permute(260*m.wQ4, [2 1 3]))
% isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
% permute((1e+3)*m.fxQ4,[2 1 3]), -fth)
% isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
% permute((1e+3)*m.fxQ4,[2 1 3]), fth)
% clim([-1.5 1.5])
% colormap jet
% colorbar
% lightangle(-45,-90)
% axis equal
% % xlim([-0.5 0.5])
% % ylim([-0.5 0.5])
% % zlim([0 0.4])
% view(45,45)
% xlabel('z')
% ylabel('x')
% zlabel('y')
% % 
% %  f2=sprintf("iso_lse_woy_w_j_%03d.fig",jcond)
% %  saveas(h1,f2)
