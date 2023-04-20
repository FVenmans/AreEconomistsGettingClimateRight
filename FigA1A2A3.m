%Make Figures A1 A2 and A3 (figure 1 for pulse of 1 GtC, 1000GtC or, 100Gt but with background RCP4.5 instead of constant background concentration)

%IN CASE OF CONSTANT BACKGROUND CONCENTRATION (Fig A1 and A2)
    % Run file FIG1_3_4, from line 1 and until GRAPHS but change Pulse=...
    %Adapt PAGE and FUND trajectory in this file below line 46-49 and name graph
%IN CASE OF BACKGROUND RCP 4.5 (Fig A3)
    %Rune file Fig1_3_4, until line 34 (but change Pulse=...)
    Epath=Epath(4,:); %RCP4.5 
    Epath=[ones(1,5)*Epath(1), Epath, ones(1,110)*Epath(86)]; %From yearly 2015 to 2100 to yearly 2010 to 2210, post 2100 at level of 2100
    % Run file Fig1_3_4, from line 40 and until GRAPHS
    %Adapt PAGE and FUND trajectory in this file below line 46-49 and name graph
%IMPORT PAGE AND FUND
    T_PAGE_1Gt          =xlsread('FUND and PAGE output for Frank 230820',1,'B4:GT4','basic');
    T_PAGE_1Gt=interp1(1:193,T_PAGE_1Gt(1,1:193), year,'spline','extrap'); 
    T_PAGE_1000Gt       =xlsread('FUND and PAGE output for Frank 230820',1,'B6:GT6','basic');
    T_PAGE_1000Gt=interp1(1:193,T_PAGE_1000Gt(1,1:193), year,'spline','extrap'); 
    T_PAGE_100Gt_RCP45  =xlsread('FUND and PAGE output for Frank 230820',1,'B8:GT8','basic');
    T_PAGE_100Gt_RCP45=interp1(1:193,T_PAGE_100Gt_RCP45(1,1:193), year,'spline','extrap');
    T_FUND_1Gt          =xlsread('FUND and PAGE output for Frank 230820',1,'B12:GT12','basic');
    T_FUND_1Gt=interp1(1:193,T_FUND_1Gt(1,1:193), year,'spline','extrap');
    T_FUND_1Gt=interp1(1:193,T_FUND_1Gt(1,1:193), year,'spline','extrap');% I do not understand why i have to extrapolate twice
    T_FUND_1000Gt       =xlsread('FUND and PAGE output for Frank 230820',1,'B14:GT14','basic');
    T_FUND_1000Gt=interp1(1:193,T_FUND_1000Gt(1,1:193), year,'spline','extrap');
    T_FUND_100Gt_RCP45  =xlsread('FUND and PAGE output for Frank 230820',1,'B16:GT16','basic'); 
    T_FUND_100Gt_RCP45=interp1(1:193,T_FUND_100Gt_RCP45(1,1:193), year,'spline','extrap');
%Fig A1 or A2 or A3 
%IRF DICE PAGE FUND Lemoine&Rudik Golosov Gerlagh&Liski
subplot(1,1,1)
    ar=area(transpose(year),transpose(T_deciles(9,:)), 'FaceColor' , [0.8 1 0.8]);
    hold on
    area(transpose(year),transpose(T_deciles(8,:)), 'FaceColor' , [0.5 0.9 0.5])
    area(transpose(year),transpose(T_deciles(7,:)), 'FaceColor' , [0 0.8 0.4])
    area(transpose(year),transpose(T_deciles(6,:)), 'FaceColor' , [0 0.6 0.3])
    area(transpose(year),transpose(T_deciles(4,:)), 'FaceColor' , [0 0.8 0.4])
    area(transpose(year),transpose(T_deciles(3,:)), 'FaceColor' , [0.5 0.9 0.5])
    area(transpose(year),transpose(T_deciles(2,:)), 'FaceColor' , [0.8 1 0.8])
    area(transpose(year),transpose(T_deciles(1,:)), 'FaceColor' , [1 1 1])
    p1=plot(year,T_JGbestfit_diff,'r','linewidth',2);
    p2=plot(year, T_FAIR_diff(1,:),'-b','linewidth',2);
    %p3=plot(year,ones(1,201)*0.00168*Pulse/3.664,'--k','linewidth',1); % cumulative emissions
    %p3=plot(year, T_FAIR_diff(4,:),'-.b','linewidth',2);        
        plot(year, T_DICE_diff(1,:),'k','linewidth',2);
    p4=scatter(year(1:10:end), T_DICE_diff(1,1:10:end),'s','MarkerEdgeColor',[0.25 0.6 1],'linewidth',1.5);
        plot(year, T_DICE2013_diff(1,:),'k','linewidth',1.5);
    p5=scatter(year(1,1:10:end), T_DICE2013_diff(1,1:10:end),'p','MarkerEdgeColor','c','linewidth',1.5);
         plot(year, T_FUND_1Gt,'k','linewidth',1.5);
    p6=scatter(year(1:10:end), T_FUND_1Gt(1,1:10:end),'v','MarkerEdgeColor','m','linewidth',1.5); 
         plot(year, T_PAGE_1Gt,'k','linewidth',1.5); 
    p7=scatter(year(1:10:end), T_PAGE_1Gt(1:10:end),'d','MarkerEdgeColor','r','linewidth',1.5);  
         plot(year, T_Golosov, 'k','linewidth',1.5);
    p8=scatter(year(1:10:end), T_Golosov(1:10:end), '+','MarkerEdgeColor',[0.9 0.55 0],'linewidth',1.5);
        plot(year, T_GL_diff, 'k', 'linewidth', 1.5);        
    p9=scatter(year(1:10:end), T_GL_diff(1:10:end), 'x','MarkerEdgeColor',[1 0.9 0.1],'linewidth', 2);
        plot(year, T_LR_diff(1,:),'k','linewidth',1.5);
    p10=scatter(year(1:10:end), T_LR_diff(1,1:10:end),'^','MarkerEdgeColor','b','linewidth',1.5);
    hold off
    legend([ar p1 p2 p4 p5 p6 p7 p8 p9 p10],'Deciles CMIP5 combinations','Best fit CMIP5 ensemble','FAIR',...
        'DICE16','DICE13','FUND','PAGE','GHKT14','GL18','LR17',...
        'Location','SouthEast','FontSize',11,'NumColumns',2) %,'Orientation','horizontal' or 'Southoutside'
    axis([0 200 0 0.0032*Pulse/3.664])  %0.32 for 100GtC pulse
    ax = gca;
    ax.YRuler.Exponent = 0;
    xlabel('years')
    ylabel('Temperature increase (°C)')
    pbaspect([1.5 1 1])
    set(gca,'fontsize',11)
    set(gcf,'PaperOrientation','landscape');
    print('-fillpage','Fig1_1GtConstantBackground', '-dpdf'); % '-fillpage' gives a too high figure