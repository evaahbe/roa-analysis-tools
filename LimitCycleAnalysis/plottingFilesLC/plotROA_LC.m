function plotROA_LC(filenames,fignum)

    cd(filenames.resultsdirectory)
    k = 0;
    file_list = dir('*.mat'); %load the file names into a list
    for j = 1:length(file_list)
        cur_file_name = file_list(j).name;
        plot_flag = 1;
        for i=1:length(filenames.filenamestoplot)
            if isempty(strfind(cur_file_name,filenames.filenamestoplot{i}))
               plot_flag = 0;
               break
            end
        end
        if plot_flag
            k = k+1;
            filestoplot_list{k} = cur_file_name;
        end
    end

    cd ..
    
    load(strcat(filenames.resultsdirectory,filestoplot_list{1}),'gamma','Qvals','numsets','sys') %just for setting the comps up

    computeROAvol = 1;

    tau_array = linspace(0,sys.Tperiod,sys.Nhp);
    tau_array_fine = linspace(0,sys.Tperiod,sys.Nhp*20);
    
    for i=1:length(tau_array_fine)
        xorbi_fine(:,i) = ppval(sys.xorb,tau_array_fine(i));    
    end


    width = 40;
    height = 28;
%     roaplot = figure('Units','centimeters',...
%     'Position',[1 0.7 width height],...
%     'PaperPositionMode','auto');
    figure(fignum)
    hold on
    colors = getCustomColorsLC();

    if sys.xdim == 2
        plot(xorbi_fine(1,:),xorbi_fine(2,:),'LineWidth',2,'Color',colors.orbit)
        legendentries{1} ='$\omega^*(\tau)$';
        
        
        domlim = max(max(abs(xorbi_fine.*[4;4]))); %this is just for a first guess
        domain = [-domlim domlim];

        Nx = 10000;

        xg = linspace(domain(1),domain(2),Nx);

        for jj = 1:length(filestoplot_list)
            disp(filestoplot_list{jj})
            
            filetoplotname = strcat(filenames.resultsdirectory,filestoplot_list{jj});
            load(filetoplotname)
            legendentries{jj+1} =  strcat(strcat(filetoplotname(strfind(filetoplotname,sys.method):end-4)));
                             

            for i = 1:sys.Nhp-1
                
                [Z_array,~] = projectOperator(sys,tau_array(i));
                xorbi = ppval(sys.xorb,tau_array(i));

                
                if length(gamma(:))~=1
                    gammai = gamma(i);
                else
                    gammai = gamma;
                end

                Qmat =reshape(Qvals(i,:),sqrt(length(Qvals(i,:))),sqrt(length(Qvals(i,:))));
                
                lypgrid = zeros(1,length(xg(:)));
                if numsets.degs.V_dU == 2
                    for j=1:length(xg(:))
                        rho1 = xg(j);
                        lypgrid(j) = [rho1]'* Qmat* [rho1]-gammai;
                    end
                elseif numsets.degs.V_dU == 4
                    for j=1:length(xg(:))
                        rho1 = xg(j);
                        lypgrid(j) = [rho1;rho1^2]'* Qmat* [rho1;rho1^2]-gammai;
                    end
                elseif numsets.degs.V_dU == 6
                    for j=1:length(xg(:))
                        rho1 = xg(j);
                        lypgrid(j) = [rho1;rho1^2;rho1^3]'* Qmat* [rho1;rho1^2;rho1^3]-gammai;
                    end

                end

                lypspline = spline(xg',lypgrid');
                
                prod_zeros = fnzeros(lypspline);

                
                xsub1 = Z_array' * (prod_zeros(1,1))+ xorbi;
                xsub2 = Z_array' * (prod_zeros(1,2))+ xorbi;
                

                if i == 1
                    plot([xsub1(1),xsub2(1)],[xsub1(2),xsub2(2)],'-', 'LineWidth',1.2,'Color',colors.roaarray{jj})
                else
                    plot([xsub1(1),xsub2(1)],[xsub1(2),xsub2(2)],'-', 'LineWidth',1.2,'Color',colors.roaarray{jj},'HandleVisibility','off')
                end  

       
                if computeROAvol == 1
                    roavol(i) = norm(xsub1-xsub2);
                    if i == sys.Nhp-1
                    	save(filetoplotname,'-append','roavol')
                    end
                end

            end
        end



    elseif sys.xdim == 3
        plot3(xorbi_fine(2,:),xorbi_fine(1,:),10*xorbi_fine(3,:),'LineWidth',2,'Color',colors.orbit)
        legendentries{1} ='$\omega^*(\tau)$';

        %title

        domlim = max(max(abs(xorbi_fine.*repmat([1;1;10],1,sys.Nhp*20)))); %this is just for a first guess
        domain = [-domlim domlim -domlim domlim]*1;

        Nx = 85;
        Ny = 85;

        xgp = linspace(domain(1),domain(2),Nx);
        ygp = linspace(domain(3),domain(4),Ny);

        [xg,yg] = meshgrid(xgp,ygp);


        contminpoints = 80;
        minpoints = 20;

        for jj = 1:length(filestoplot_list)

            filetoplotname = strcat(filenames.resultsdirectory,filestoplot_list{jj});
            load(filetoplotname)
            legendentries{jj+1} =  strcat(strcat(filetoplotname(strfind(filetoplotname,sys.method):end-4)));
            
%             for i=1:length(tau_array)
%                 Z_array(:,:,i) = projectOperator(sys,tau_array(i));    
%             end

            for i = 1:sys.Nhp-1
                xorbi = ppval(sys.xorb,tau_array(i));
                Z_array = projectOperator(sys,tau_array(i));


                Qmat =reshape(Qvals(i,:),sqrt(length(Qvals(i,:))),sqrt(length(Qvals(i,:))));

                lypgrid = zeros(1,length(xg(:)));
                if numsets.degs.V_dU == 2
                    for j=1:length(xg(:))
                        rho1 = xg(j);
                        rho2 = yg(j);
                        lypgrid(j) = [rho1;rho2]'* Qmat* [rho1;rho2];
                    end
                elseif numsets.degs.V_dU == 4
                    for j=1:length(xg(:))
                        rho1 = xg(j);
                        rho2 = yg(j);
                        lypgrid(j) = [rho1;rho2;rho1^2;rho1*rho2;rho2^2]'* Qmat* [rho1;rho2;rho1^2;rho1*rho2;rho2^2];
                    end
                elseif numsets.degs.V_dU == 6
                    for j=1:length(xg(:))
                        rho1 = xg(j);
                        rho2 = yg(j);
                        lypgrid(j) = [rho1;rho2;rho1^2;rho1*rho2;rho2^2;rho1^3;rho1^2*rho2;rho1*rho2^2;rho2^3]'* Qmat* [rho1;rho2;rho1^2;rho1*rho2;rho2^2;rho1^3;rho1^2*rho2;rho1*rho2^2;rho2^3];
                    end

                end

                if length(gamma(:))~=1
                    gammai = gamma(i);
                else
                    gammai = gamma;
                end

                lypgrid = reshape(lypgrid,size(xg));
                vlyp = [gammai, gammai];
                ctmatrixcalc = contourc(xgp,ygp,lypgrid,vlyp);  
                ctmatrixcalc = ctmatrixcalc(:,2:end);
                ctmatrixcalc = (abs(ctmatrixcalc)<=domlim).*ctmatrixcalc + ((ctmatrixcalc)>domlim).*domlim  -((ctmatrixcalc)<-domlim).*domlim;
                
                xsub = zeros(sys.xdim,length(ctmatrixcalc(1,:))); 
                for j = 1:length(ctmatrixcalc(1,:))
                    xsub(:,j) = Z_array' * ctmatrixcalc(:,j);
                end

                subcoords = xsub + xorbi;
            
                ctmatrixcalc_clean = ctmatrixcalc;

                
                if i == 1
                    plot3(subcoords(2,:),subcoords(1,:),10*subcoords(3,:),'-', 'LineWidth',1.2,'Color',colors.roaarray{jj})
                else
                    plot3(subcoords(2,:),subcoords(1,:),10*subcoords(3,:),'-', 'LineWidth',1.2,'Color',colors.roaarray{jj},'HandleVisibility','off')
                end             
                
                contcoords =[subcoords(1,:),subcoords(2,:),10*subcoords(3,:)];
                contourcoords{i} = contcoords;
                
                
                
                if computeROAvol == 1

                    continterv = linspace(0,1,length(ctmatrixcalc_clean(1,:)));
                    contspline = spline(continterv,ctmatrixcalc_clean);
                    contnormint = linspace(0,1,contminpoints);
                    for j = 1:contminpoints
                        contline(:,j) = ppval(contspline,contnormint(j));
                    end
                    roavol(i) = polyarea(contline(2,:),contline(1,:));
                    if i == sys.Nhp-1
                        save(filetoplotname,'-append','roavol')
                        save(filetoplotname,'-append','contourcoords')
                    end
                end
                

            end
        end

     else
         error('Too many dimensions to plot')
     end


    legend(legendentries,...
        'interpreter','latex',...
        'FontSize',20,...
        'Location','NorthEast')

    set(gca,...
    'LineWidth', 1.5,...             
    'FontSize',20)
  
    if strcmp(sys.system, 'Kite')
         ylabel({'Elevation $\mathbf{\theta}$ [rad]'},...
        'interpreter','latex',...
        'FontWeight','bold',...
        'FontSize',35)

        xlabel('Azimuth $\mathbf{\phi}$ [rad]',...
        'interpreter','latex',...
        'FontWeight','bold',...
        'FontSize',35)

        zlabel('Orientation $\mathbf{\gamma}$ [rad]',...
        'interpreter','latex',...
        'FontWeight','bold',...
        'FontSize',35)
        %axis([-0.8 0.8 0.1 0.6 -3 3]);
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


