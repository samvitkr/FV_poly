close all
clear
jcond=156;
lt=0.05;
%fth=0.1;
%ft=sprintf("velfield_lse_voz_j_%d.mat",jcond)
%ft=sprintf("velfield_lse_vwoy_j_%d.mat",jcond)
ft=sprintf('../data/velfield_lse_ddfilter_j_%03d.mat',jcond);
m=matfile(ft)
ftu=sprintf('../data/velfield_lse_uufilter_j_%03d.mat',jcond);
mu=matfile(ftu)
x1=150;
y1=150;
x2=1000;
y2=350;
h1=figure('OuterPosition',...
    [x1 y1 x2 y2]);
subplot(1,2,1)
isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
permute(m.Q2Q,[2 1 3]), lt)
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
isosurface( permute(mu.Z,[2 1 3]), permute(mu.X,[2 1 3]), permute(mu.Y,[2 1 3]),...
permute(mu.Q2Q,[2 1 3]), lt)
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

%subplot(1,3,3)
%isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
%permute(m.Q2lambda2,[2 1 3]), lt, permute(260*m.wQ2, [2 1 3]))
% isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
% permute((1e+3)*m.fzQ2,[2 1 3]), -fth)
% isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
% permute((1e+3)*m.fzQ2,[2 1 3]), fth)
%clim([-1 1])
%colormap jet
%colorbar
%lightangle(-45,-90)
%axis equal
%ylim([-1 1])
%xlim([-0.5 0.5])
%% zlim([0 0.4])
%view(45,45)
%xlabel('z')
%ylabel('x')
%zlabel('y')
%grid on
f1=sprintf("iso_lse_vel_Q_filter_j_%03d.fig",jcond)
saveas(h1,f1)
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
