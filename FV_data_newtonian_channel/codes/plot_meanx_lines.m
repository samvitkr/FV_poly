clear
close all
jcond=205;
jc=jcond-110;
lx=4*pi;
nx=512;
x=lx*[0:nx-1]./nx-lx/2;
mm=matfile('../data/mean_profiles.mat');
ut=mm.ret/mm.re;
load('../data/ygrid.mat')
y=yCheb+1;
yc=y(jcond);
y=y(111:end);
%
xl1=-0.4;
xl2=0.401;
yl1=0;
yl2=0.81;
%snp=sprintf("../data/velgrad_v_profile_j_%03d.mat",jcond);
%m=matfile(snp);
% sd=sprintf("xslices_j_%03d_down_2",jcond);
% su=sprintf("xslices_j_%03d_up_2",jcond);
% vdw = VideoWriter(sd,'MPEG-4');
% vdw.FrameRate=5;
% open(vdw)
% vuw = VideoWriter(su,'MPEG-4');
% vuw.FrameRate=5;
% open(vuw)
    x1=100;
    y1=40;
    x2=900;
    y2=700;
    %%
    fd=figure('OuterPosition',...
        [x1 y1 x2 y2]);
for i =210:280
    
    islice=i;
    xi=x(islice);
    fn=sprintf("../data/xslices/velgrad_xslice_lsevpn_i_%03d_j_%03d.mat",islice,jcond)
    load(fn)
    td=sprintf("x=%.02f",xi)
    % picd=sprintf("plotd_i_%03d_j_%03d.png",islice,jcond)
    % picu=sprintf("plotu_i_%03d_j_%03d.png",islice,jcond)

    %%

    %
    subplot(3,3,1)
    pcolor(Z,Y,vd./ut)
    shading flat
    axis equal
    xlim([xl1 xl2])
    ylim([yl1 yl2])
    clim([-1e-1 1e-1])
    title('(a) $v/u_{\tau}$' ,'Interpreter','latex')
    ylabel('y')
    xlabel('z')
    xline(0)
    yline(yc)
    colorbar

    subplot(3,3,2)
    pcolor(Z,Y,wd./ut)
    shading flat
    axis equal
    xlim([xl1 xl2])
    ylim([yl1 yl2])
    title('(b) $w/u_{\tau}$' ,'Interpreter','latex')
    ylabel('y')
    xlabel('z')
    xline(0)
    yline(yc)
    colorbar
    clim([-1e-1 1e-1])
    colormap redblue

    %
    [startX,startY]=meshgrid(0.4,0.05:0.05:0.8);
    %[startX2,startY2]=meshgrid(-0.27,0.05:0.02:0.2);

    subplot(3,3,3)
    vertsv = stream2(Z,Y,wd,vd,-startX,startY);
    vertsv2 = stream2(Z,Y,wd,vd,startX,startY);
    hold on
    pcolor(Z,Y,ud./ut)
    l1=streamline(vertsv);
    l2=streamline(vertsv2);
    shading flat
    hold off
    axis equal
    xlim([xl1 xl2])
    ylim([yl1 yl2])
    title("$(c)u'/u_{\tau}$"  ,'Interpreter','latex')
    ylabel('y')
    xlabel('z')
    xline(0)
    yline(yc)
    colorbar
    clim([-0.2 0.2])    
set(l1, 'Color', 'k'); 
set(l2, 'Color', 'k'); 

        subplot(3,3,4)
    pcolor(Z,Y,-ozd./ut)
    shading flat
    axis equal
    xlim([xl1 xl2])
    ylim([yl1 yl2])
    title('(d) $-\omega_z/(u_{\tau}/H)$' ,'Interpreter','latex')
    ylabel('y')
    xlabel('z')
    xline(0)
    yline(yc)
    colorbar
    clim([-250 250])

    colormap redblue
    subplot(3,3,5)
    pcolor(Z,Y,oyd./ut)
    shading flat
    axis equal
    xlim([xl1 xl2])
    ylim([yl1 yl2])
    clim([-5 5])
    title('(e) $\omega_y/(u_{\tau}/H)$' ,'Interpreter','latex')
    ylabel('y')
    xlabel('z')
    xline(0)
    yline(yc)
    colorbar


    %
    [startX,startY]=meshgrid(0.4,0.02:0.02:0.8);
    [startX2,startY2]=meshgrid(-0.27,0.05:0.02:0.2);

    subplot(3,3,6)
    verts = stream2(Z,Y,ozd,oyd,startX,startY);
    verts2 = stream2(Z,Y,-ozd,-oyd,startX2,startY2);
    hold on
     pcolor(Z,Y,oxd./ut)
    shading flat
    l1=streamline(verts);
    %l2=streamline(verts2);
    hold off
    axis equal
    xlim([xl1 xl2])
    ylim([yl1 yl2])
    title('(f) $\omega_x/(u_{\tau}/H)$' ,'Interpreter','latex')
    ylabel('y')
    xlabel('z')
    xline(0)
    yline(yc)
    colorbar
        clim([-5 5])

    sgtitle(td)
set(l1, 'Color', 'k'); 


    subplot(3,3,7)
   % pcolor(Z,Y,-(ozd.*vd)./ut^2)
   pcolor(Z,Y,-(vozdc)./ut^2)
   shading flat
    axis equal
    xlim([xl1 xl2])
    ylim([yl1 yl2])
    clim([-2 2])
    title('(g) $v \omega_z /(-u_{\tau}^2/H)$' ,'Interpreter','latex')
    ylabel('y')
    xlabel('z')
    xline(0)
    yline(yc)
    colorbar

    subplot(3,3,8)
    %pcolor(Z,Y,(oyd.*wd)./ut^2)
    pcolor(Z,Y,(woydc)./ut^2)    
    shading flat
    axis equal
    xlim([xl1 xl2])
    ylim([yl1 yl2])
    title('(h) $-w\omega_y/(-u_{\tau}^2/H)$' ,'Interpreter','latex')
    ylabel('y')
    xlabel('z')
    xline(0)
    yline(yc)
    colorbar
    clim([-2 2])
    colormap redblue

    subplot(3,3,9)
    pcolor(Z,Y,ld)
    shading flat
    axis equal
    xlim([xl1 xl2])
    ylim([yl1 yl2])
   clim([-5e-3 5e-3])
    title('(e) $\lambda_2/(U^2/H^2)$' ,'Interpreter','latex')
    ylabel('y')
    xlabel('z')
    xline(0)
    yline(yc)
    colorbar

    %exportgraphics(fd,picd,'BackgroundColor','none')
    frame = getframe(fd);
    writeVideo(vdw,frame)
    clf
    %%

   % 
   % fu=figure('OuterPosition',...
   %     [x1 y1 x2 y2]);

    subplot(3,3,1)
    pcolor(Z,Y,vu./ut)
    shading flat
    axis equal
    xlim([xl1 xl2])
    ylim([yl1 yl2])
    clim([-1e-1 1e-1])
    title('(a) $v/u_{\tau}$' ,'Interpreter','latex')
    ylabel('y')
    xlabel('z')
    xline(0)
    yline(yc)
    colorbar

    subplot(3,3,2)
    pcolor(Z,Y,wu./ut)
    shading flat
    axis equal
    xlim([xl1 xl2])
    ylim([yl1 yl2])
    title('(b) $w/u_{\tau}$' ,'Interpreter','latex')
    ylabel('y')
    xlabel('z')
    xline(0)
    yline(yc)
    colorbar
    clim([-0.1 0.1])
    colormap redblue

    %
    [startX,startY]=meshgrid(0.4,0.05:0.05:0.8);
    %[startX2,startY2]=meshgrid(-0.27,0.05:0.02:0.2);



    subplot(3,3,3)
    vertsv = stream2(Z,Y,wu,vu,-startX,startY);
    vertsv2 = stream2(Z,Y,wu,vu,startX,startY);
    hold on
    pcolor(Z,Y,uu./ut)
    l1=streamline(vertsv);
    l2=streamline(vertsv2);
    shading flat
    hold off
    axis equal
    xlim([xl1 xl2])
    ylim([yl1 yl2])
    title("$(c)u'/u_{\tau}$"  ,'Interpreter','latex')
    ylabel('y')
    xlabel('z')
    xline(0)
    yline(yc)
    colorbar
    clim([-0.2 0.2])
set(l1, 'Color', 'k'); 
set(l2, 'Color', 'k'); 

subplot(3,3,4)
    pcolor(Z,Y,-ozu./ut)
    shading flat
    axis equal
    xlim([xl1 xl2])
    ylim([yl1 yl2])
    title('(e) $-\omega_z/(u_{\tau}/H)$' ,'Interpreter','latex')
    ylabel('y')
    xlabel('z')
    xline(0)
    yline(yc)
    colorbar
    clim([-250 250])

    colormap redblue

    subplot(3,3,5)
    pcolor(Z,Y,oyu./ut)
    shading flat
    axis equal
    xlim([xl1 xl2])
    ylim([yl1 yl2])
    clim([-5 5])
    title('(d) $\omega_y/(u_{\tau}/H)$' ,'Interpreter','latex')
    ylabel('y')
    xlabel('z')
    xline(0)
    yline(yc)
    colorbar


    %
    [startX,startY]=meshgrid(0.4,0.02:0.02:0.8);
    [startX2,startY2]=meshgrid(-0.27,0.05:0.02:0.2);

    subplot(3,3,6)
    verts = stream2(Z,Y,ozu,oyu,startX,startY);
    %verts2 = stream2(Z,Y,-ozu,-oyu,startX2,startY2)
    hold on
    pcolor(Z,Y,oxu./ut)
    l1=streamline(verts);
    %l2=streamline(verts2);
    shading flat
    hold off
    axis equal
    xlim([xl1 xl2])
    ylim([yl1 yl2])
    title('$(f)\omega_x/(u_{\tau}/H)$' ,'Interpreter','latex')
    ylabel('y')
    xlabel('z')
    xline(0)
    yline(yc)
    colorbar
    clim([-5 5])
    sgtitle(td)
set(l1, 'Color', 'k'); 

    subplot(3,3,7)
    %pcolor(Z,Y,-(ozu.*vu)./ut^2)
    pcolor(Z,Y,-(vozuc)./ut^2)
    shading flat
    axis equal
    xlim([xl1 xl2])
    ylim([yl1 yl2])
    clim([-2 2])
    title('(g) $v \omega_z /(-u_{\tau}^2/H)$' ,'Interpreter','latex')
    ylabel('y')
    xlabel('z')
    xline(0)
    yline(yc)
    colorbar

    subplot(3,3,8)
  %  pcolor(Z,Y,(oyu.*wu)./ut)
    pcolor(Z,Y,(woyuc)./ut^2)
    shading flat
    axis equal
    xlim([xl1 xl2])
    ylim([yl1 yl2])
    title('(h) $-w\omega_y/(-u_{\tau}^2/H)$' ,'Interpreter','latex')
    ylabel('y')
    xlabel('z')
    xline(0)
    yline(yc)
    colorbar
    clim([-2 2])
    colormap redblue

     subplot(3,3,9)
    pcolor(Z,Y,lu)
    shading flat
    axis equal
    xlim([xl1 xl2])
    ylim([yl1 yl2])
   clim([-5e-3 5e-3])
    title('(e) $\lambda_2/(U^2/H^2)$' ,'Interpreter','latex')
    ylabel('y')
    xlabel('z')
    xline(0)
    yline(yc)
    colorbar
%     %exportgraphics(fu,picu,'BackgroundColor','none')
     frameu = getframe(fd);
    writeVideo(vuw,frameu)
    clf
end
close(vdw)
close(vuw)