function colors = getCustomColorsLC()

    colors = struct();
    colors.mycolors = [  
            0.8500,    0.3250,    0.0980;...
            0,         0.4470,    0.7410;...
            0.9290,    0.6940,    0.1250;...
            0.4660,    0.6740,    0.1880;...
            0.4940,    0.1840,    0.5560;...
            0.3010,    0.7450,    0.9330;...
            0.6350,    0.0780,    0.1840;...
            0.8,       0.8,       0.6;...
            0.2,       0.9,       0.6;...
            0,         0.1,       0.1];
    
    colors.orbit = [0 0 0];
    colors.roa1 = [0 0.4 1];
    colors.roa2 = [0.4 0.4 0];
    colors.traj1 = [1 0.2 0];

    colors.roaarray{1} = colors.mycolors(1,:);
    colors.roaarray{2} = colors.mycolors(2,:);
    colors.roaarray{3} = colors.mycolors(3,:);
    colors.roaarray{4} = colors.mycolors(4,:);
    colors.roaarray{5} = colors.mycolors(5,:);
    colors.roaarray{6} = colors.mycolors(6,:);
    colors.roaarray{7} = colors.mycolors(7,:);
    colors.roaarray{8} = colors.mycolors(8,:);
    colors.roaarray{9} = colors.mycolors(9,:);
    colors.roaarray{10} = colors.mycolors(10,:);
    
end
