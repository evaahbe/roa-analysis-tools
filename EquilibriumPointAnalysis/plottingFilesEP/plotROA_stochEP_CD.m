% this script only works for 2D stochastic systems!

function plotROA_stochEP_CD(filenames)

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
    
    width = 28;
    height = 28;
    roaplot = figure('Units','centimeters',...
    'Position',[1 0.7 width height],...
    'PaperPositionMode','auto');
    hold on
    for jj = 1:length(filestoplot_list)
        
        filetoplotname = strcat(filenames.resultsdirectory,filestoplot_list{jj});
        load(filetoplotname)
        
        colors = getCustomColors();

        domain = [ -2*sqrt(1/Q0(1,1)) 2*sqrt(1/Q0(1,1)) -2*sqrt(1/Q0(2,2)) 2*sqrt(1/Q0(2,2))]; %this is a first rough estimate, adjust if needed.
        Nx = 200;
        Ny = 200;
        xgp = linspace(domain(1),domain(2),Nx);
        ygp = linspace(domain(3),domain(4),Ny);
        [xg,yg] = meshgrid(xgp,ygp);
        
        lypgrid = zeros(1,length(xg(:)));
        if numsetsRE.degs.V_dU == 2
            for j=1:length(xg(:))
                x01 = xg(j);
                x02 = yg(j);
                lypgrid(j) = [x01;x02]'* Q0* [x01;x02];
            end
        elseif numsetsRE.degs.V_dU == 4
            for j=1:length(xg(:))
                x01 = xg(j);
                x02 = yg(j);
                lypgrid(j) = [x01;x02;x01^2;x01*x02;x02^2]'* Q0* [x01;x02;x01^2;x01*x02;x02^2];
            end
        elseif numsetsRE.degs.V_dU == 6
            for j=1:length(xg(:))
                x01 = xg(j);
                x02 = yg(j);
                lypgrid(j) = [x01;x02;x01^2;x01*x02;x02^2;x01^3;x01^2*x02;x01*x02^2;x02^3]'* Q0* [x01;x02;x01^2;x01*x02;x02^2;x01^3;x01^2*x02;x01*x02^2;x02^3];
            end
        else 
            error('Dimension not implemented.')
        end
        lypgrid = reshape(lypgrid,size(xg));
        vlyp = [1 1]; %if alpha=1, otherwise insert alpha value
        contourmatrix = contourc(xgp,ygp,lypgrid,vlyp);
        contourmatrix = contourmatrix(:,2:end);
        domlim = max(abs(domain));
        contourmatrix = (abs(contourmatrix)<=domlim).*contourmatrix + ((contourmatrix)>domlim).*domlim  -((contourmatrix)<-domlim).*domlim;
            
        plot(contourmatrix(1,:)+xiep0(1),contourmatrix(2,:)+xiep0(1),'LineWidth',1.,'Color',colors.mycolors(jj,:))
        %legend entries, only showing diagonal elements
        plot_name = strcat(strcat(filetoplotname(strfind(filetoplotname,'ults_')+5:strfind(filetoplotname,'ults_')+7)));
        h1_name = strcat(strcat(filetoplotname(strfind(filetoplotname,'h1'):strfind(filetoplotname,'h1')+1)));
        h2_name = strcat(strcat(filetoplotname(strfind(filetoplotname,'h2'):strfind(filetoplotname,'h2')+1)));
        legendentries{jj} =  strcat('$\partial(V) =$ ',num2str(numsetsRE.degs.V_dU),', ',plot_name,', ',h1_name,h2_name);
        
        % OPTIONAL: save contour data for Monte Carlo testing
        save(strcat(filenames.resultsdirectory,filestoplot_list{jj}), '-append','contourmatrix');
    end   
    
    
    plot(xiep0(1), xiep0(2), 'ko','LineWidth',1.)
    legendentries{jj+1} = '$x_{EP}$';

    set(gca,...
    'LineWidth', 1.2,...
    'Units','normalized',...                  
    'FontSize',20)
    title('ROA estimates','FontSize',20)
    legend(legendentries,'LineWidth',1.2,'FontSize',23,'NumColumns',1,'interpreter','latex')
    xlabel('$x_1$','FontSize',35,'interpreter','latex');
    ylabel('$x_2$','FontSize',35,'interpreter','latex');  

end




    
    

    
    
