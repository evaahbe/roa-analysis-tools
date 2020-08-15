function plotHyperplanes_LC(sys)

    tau_array = linspace(0,sys.Tperiod,sys.Nhp);

    for i=1:length(tau_array)
        xorbi(:,i) = ppval(sys.xorb,tau_array(i));    
    end
    
    tau_array_fine = linspace(0,sys.Tperiod,sys.Nhp*20);

    for i=1:length(tau_array_fine)
        xorbi_fine(:,i) = ppval(sys.xorb,tau_array_fine(i));    
    end

    width = 40;
    height = 28;
%     roaplot = figure('Units','centimeters',...
%     'Position',[1 0.7 width height],...
%     'PaperPositionMode','auto');
figure(2)
    hold on
    colors = getCustomColorsLC();
    

    if sys.xdim == 2
        plot(xorbi_fine(1,:),xorbi_fine(2,:),'LineWidth',2,'Color',colors.orbit)
        legendentries{1} ='$\omega^*(\tau)$';
              

        for i = 1:sys.Nhp-1

            Z_array(:,:,i) = projectOperator(sys,tau_array(i));  
            
            rhoval = 2;

            xsub(:,1) = Z_array(:,:,i)' * (+rhoval);
            if sys.trafo ==2
                xsub(:,2) = -xorbi(:,i)+sys.xc;
            else
                xsub(:,2) = Z_array(:,:,i)' * (-rhoval);
            end

            subcoords = xsub + xorbi(:,i);
            if i == 1
                plot(subcoords(1,:),subcoords(2,:),'-', 'LineWidth',1.,'Color',colors.roaarray{1})
            else
                plot(subcoords(1,:),subcoords(2,:),'-', 'LineWidth',1.,'Color',colors.roaarray{1},'HandleVisibility','off')
                
            end
            
        end
        
        
        legendentries{2} ='$S(\tau)$';
        if sys.trafo ==2
            plot(sys.xc(1),sys.xc(2),'ko','LineWidth',1.)
            legendentries{3} ='$x_c$';
        end


    elseif sys.xdim == 3
        if strcmp(sys.system, 'Kite')
            plot3(xorbi_fine(1,:),xorbi_fine(2,:),10*xorbi_fine(3,:),'LineWidth',2,'Color',colors.orbit)
            legendentries{1} ='$\omega^*(\tau)$';

            %title

            spec_taulist = linspace(0.5,sys.Tperiod,7);
            domlim = max(max(abs(xorbi.*repmat([1;1;1],1,sys.Nhp))))*1.5; %this is just for a first guess
            domain = [-domlim domlim -domlim domlim];
        else
            plot3(xorbi(1,:),xorbi(2,:),xorbi(3,:),'LineWidth',2,'Color',colors.orbit)
            legendentries{1} ='$\omega^*(\tau)$';

            %title

            spec_taulist = linspace(0.1,sys.Tperiod/2,5);
            domlim = max(max(abs(xorbi.*repmat([1;1;1],1,sys.Nhp)))); %this is just for a first guess
            domain = [-domlim domlim -domlim domlim];
        end

        Nx = 5;
        Ny = 5;

        xgp = linspace(domain(1),domain(2),Nx);
        ygp = linspace(domain(3),domain(4),Ny);

        [xg,yg] = meshgrid(xgp,ygp);

        C = ones(size(xg));
        for i = 1:length(spec_taulist)-1
            tau = spec_taulist(i);
            vi = ppval(sys.v,tau);
            xref = ppval(sys.xorb,tau);
            xz = -(vi(1) * (xg - xref(1)) + vi(2) * (yg - xref(2)))/vi(3) + xref(3);
            surf(xg,yg,10*xz,C,'FaceAlpha',0.5,'EdgeAlpha',0)
            disp(max(xz))
            plot3(xref(1,:),xref(2,:),10*xref(3,:),'ko','Color',colors.orbit)
        end

     else
         error('Too many dimensions to plot')
     end


    legend(legendentries,...
        'interpreter','latex',...
        'FontSize',20,...
        'Location','NorthEast')

    set(gca,...
    'LineWidth', 1.,...
    'Units','normalized',...                  
    'FontSize',20)
  
    if strcmp(sys.system, 'Kite')
         xlabel({'Elevation $\mathbf{\theta}$ [rad]'},...
        'interpreter','latex',...
        'FontWeight','bold',...
        'FontSize',35)

        ylabel('Azimuth $\mathbf{\phi}$ [rad]',...
        'interpreter','latex',...
        'FontWeight','bold',...
        'FontSize',35)

        zlabel('Orientation $\mathbf{\gamma}$ [rad]',...
        'interpreter','latex',...
        'FontWeight','bold',...
        'FontSize',35)
        axis([ 0. 0.7 -1 1 -4 4]);
    else
        ylabel('$x_2$',...
        'interpreter','latex',...
        'FontWeight','bold',...
        'FontSize',30)

        xlabel('$x_1$',...
        'interpreter','latex',...
        'FontWeight','bold',...
        'FontSize',30)

        zlabel('$x_3$',...
        'interpreter','latex',...
        'FontWeight','bold',...
        'FontSize',30)
    end
end


