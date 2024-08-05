close all
clear
load('../data/ygrid.mat')
nx=512;
nz=384;
lx=4*pi;
lz=2*pi;
xp=lx*[0:nx-1]/nx;
zp=lz*[0:nz-1]/nz;
yp=yCheb'+1;
yp=yp(111:end);
re=4667;
[X,Z,Y]=meshgrid(xp,zp,yp);
lt=1;
t=40000;
%fth=0.1;
%ft=sprintf("velfield_lse_voz_j_%d.mat",jcond)
%ft=sprintf("velfield_lse_vwoy_j_%d.mat",jcond)
ft=sprintf('../data/lambda_ddfilter_%07d.mat',t);
m=matfile(ft)
Qm=mean(m.Q(1:end,1:end,111:end)  ,[1 2]);
Qrd=rms(m.Q(1:end,1:end,111:end)-Qm,[1,2]);

ftu=sprintf('../data/lambda_uufilter_%07d.mat',t);
mu=matfile(ftu)
Qm=mean(mu.Q(1:end,1:end,111:end)  ,[1 2]);
Qru=rms(mu.Q(1:end,1:end,111:end)-Qm,[1,2]);

x1=150;
y1=150;
x2=700;
y2=350;
h1=figure('OuterPosition',...
    [x1 y1 x2 y2]);
isosurface( permute(Z,[2 1 3]), permute(X,[2 1 3]), permute(Y,[2 1 3]),...
permute(m.Q(1:end,1:end,111:end)./Qrd,[2 1 3]), lt)
%isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
%permute((1e+3)*m.fxQ2,[2 1 3]), -fth)
%isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
%permute((1e+3)*m.fxQ2,[2 1 3]), fth)
%clim([-1 1])
colormap jet
colorbar
lightangle(-45,-90)
axis equal
%ylim([-1 1])
%xlim([-0.5 0.5])
% zlim([0 0.4])
view(45,45)
xlabel('z')
ylabel('x')
zlabel('y')
grid on

f1=sprintf("iso_Q_ddfilter_%07d.fig",t)
saveas(h1,f1)


close all

h2=figure('OuterPosition',...
    [x1 y1 x2 y2]);

isosurface( permute(Z,[2 1 3]), permute(X,[2 1 3]), permute(Y,[2 1 3]),...
permute(mu.Q(1:end,1:end,111:end)./Qru,[2 1 3]), lt)
%isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
%permute((1e+3)*m.fyQ2,[2 1 3]), -fth)
%isosurface( permute(m.Z,[2 1 3]), permute(m.X,[2 1 3]), permute(m.Y,[2 1 3]),...
%permute((1e+3)*m.fyQ2,[2 1 3]), fth)
%clim([-1 1])
colormap jet
colorbar
lightangle(-45,-90)
axis equal
%ylim([-1 1])
%xlim([-0.5 0.5])
% zlim([0 0.4])
view(45,45)
xlabel('z')
ylabel('x')
zlabel('y')
grid on

f2=sprintf("iso_Q_uufilter_%07d.fig",t)
saveas(h2,f2)
