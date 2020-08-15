function plotROAvolume_LC(filenames,fignum)

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
    
    load(strcat(filenames.resultsdirectory,filestoplot_list{1}),'sys') %just for setting the comps up

    tau_array = linspace(1,sys.Nhp-1,sys.Nhp-1);
    
    colors = getCustomColorsLC();

    figure(fignum)
    hold on

    for jj = 1:length(filestoplot_list)
            
        filetoplotname = strcat(filenames.resultsdirectory,filestoplot_list{jj});
        load(filetoplotname)
        legendentries{jj} =  strcat(strcat(filetoplotname(strfind(filetoplotname,sys.method):end-4)));
        
        roavol_resort = roavol;
            
        plot(tau_array(:), roavol_resort(:),'LineWidth',2,'Color',colors.roaarray{jj})

    end
        
        
       
    set(gca,...
    'LineWidth', 1.5,...
    'Units','normalized',...                     
    'FontWeight','bold',...
    'FontSize',20)

     ylabel({'ROA volume'},...
    'interpreter','latex',...
    'FontWeight','bold',...
    'FontSize',35)

    xlabel('$i$',...
    'interpreter','latex',...
    'FontWeight','bold',...
    'FontSize',35)


    legend(legendentries,...
        'interpreter','latex',...
    'FontSize',15,...
    'Location','NorthEast')
end
