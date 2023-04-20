%CODE FOR FIG 2 and Appendix E
%Graphs for emission scearios corresponding to RCP's
%ppm and GtC is converted to GtCO2 (above steadystate) unless indicated by the name of the matrix
close all; clear;
cd 'C:\Users\Frank\Google Drive\Research Projects\Critique DICE\GraphsDecayInertia';%Change path here
[RCP,txt_RCP,Table_RCP] = xlsread('iamc_db_Emissions.xlsx',1,'E1:O8','basic'); % does not work without 'basic'
RCP(2:8,:)=RCP(2:8,:)/1000; % convert Mt into GtCO2
[J,txt_J,Table_J] = xlsread('decay and inertia parameters.xlsx',1,'A1:H19','basic');
[G,txt_G,Table_G] = xlsread('decay and inertia parameters.xlsx',3,'A2:G20','basic');
[CONC] = xlsread('iamc_db_CO2ppm.xlsx',1,'F1:P8','basic');
[TEMP] = xlsread('iamc_db_Temp.xlsx',1,'F1:P8','basic');
[F_TOT] = xlsread('iamc_db_TotalF.xlsx',1,'F1:P8','basic');
[F_CO2] = xlsread('iamc_db_CO2F.xlsx',1,'F1:P8','basic');
[FUND_cum]=xlsread('FUND 3.11.xlsx',3,'E36:CQ43','basic');
[PAGE_cum]=xlsread('Page Cycle SD Sept 2019.xlsx',2,'D13:CP20','basic');
[IPCC_cum]=xlsread('IPCC SPM10',2,'C1:L8','basic');
    IPCC_cum(1:2:7,:)=IPCC_cum(1:2:7,:)-0.28; %IPCC graph has initial temperature 0.28°C too high, see Millar, R. J., Fuglestvedt, J. S., Friedlingstein, P., Rogelj, J., Grubb, M. J., Matthews, H. D., … Allen, M. R. (2017). Emission budgets and pathways consistent with limiting warming to 1.5 ? C. Nature Geoscience
    IPCC_cum(2:2:8,:)=(IPCC_cum(2:2:8,:)+125)*3.664; %IPCC cumulative emissions are 420GtC in IPCCAR5, but observed record is 545, see Millar, R. J., Fuglestvedt, J. S., Friedlingstein, P., Rogelj, J., Grubb, M. J., Matthews, H. D., … Allen, M. R. (2017). Emission budgets and pathways consistent with limiting warming to 1.5 ? C. Nature Geoscience.
    FUND_cum(2:2:8,:)=FUND_cum(2:2:8,:)*3.664;
    PAGE_cum(2:2:8,:)=PAGE_cum(2:2:8,:)*3.664;
F_NONCO2=F_TOT-F_CO2;
load('ALPHA_RCP.mat')
T0atm=0.85;
T0ocean=0.28; % Initial Value from FAIR
MequilibriumGtCO2=588*3.664; %background GtCO2 stock in 2015 
Climatesensitivity=3.1; %value in DICE
Ecum2015=571*3.664; %Fair excel file has571. Richard Millar in Nature has 545GtC for 2015
firstyear=2015;
lastyear=2100;
year=firstyear:1:lastyear ;
Decennia=RCP(1,:);
Decennia2005=[2005 2010:10:2100];
txt_RCP(9,1)={'Constant 2015 emissions'};RCP(9,:)=39.15;F_NONCO2(9,:)=0.01811;CONC(9,:)=missing;TEMP(9,:)=missing; 
M_DICE=zeros(8,86);M_FAIR=zeros(8,86);M_Joosbestfit=zeros(8,86);M_GL=zeros(8,86);M_Golosov=zeros(8,86);M_LR=zeros(8,86);M_GL=zeros(8,86);
T_DICE=zeros(8,86);T_FAIR=zeros(8,86);T_JoosGeofbestfit=zeros(8,86);T_Golosov=zeros(8,86);T_LR=zeros(8,86);T_GL=zeros(8,86);
X2100=zeros(8,7);
%Calculate cumulative emissions
Epath=zeros(8,86);
for i=1:8 ;Epath(i,:)=interp1(Decennia,RCP(i+1,:),year, 'spline'); end
Ecum_yearly=Ecum2015+cumsum(Epath,2);
CUM=horzcat(Ecum_yearly(:,1),Ecum_yearly(:,6),Ecum_yearly(:,16),Ecum_yearly(:,26),Ecum_yearly(:,36),Ecum_yearly(:,46),Ecum_yearly(:,56),Ecum_yearly(:,66),Ecum_yearly(:,76),Ecum_yearly(:,86));
J(19,:)=[0 1 0 0 0.0138 0 0]; %Lemoine and Rudik
J(20,:)=[0.2 0.314 0.486 0	0.00228 100	0]; %Golosov, 
J(21,:)=[0.163 0.184 0.449 0 -0.1*log(1-0.074) -0.1*log(1-0.47) 0]; %Gerlagh & Liski
%G(19,:)=[131.4 1 1.196 0 0 0]; %Lemoine and Rudik
G(19,:)=[1/0.0091/0.809 1 1/0.809 0 0 0];
epsilon=-1/10*log(1-0.183);
G(20,:)=[1.13/epsilon 1 1.13 0 0 0];%Gerlagh and Liski, using lambda from Geoffroy 
%Mi0=[0.5294;0.3433;0.1109;0.01622]*ones(1,18)*263*3.664; % proportions received from Richard Millar,  DICE=851-588=263
Mi0=[0.473184713605214;0.359104481245797;0.143649591297395;0.0240612138515941]*ones(1,18)*263*3.664;% FAIR calculated from historic emissions in InitialValuesFAIR_GL_Golosov.
Mi0(:,19)=[0 ;1; 0; 0]*263*3.6644; %Lemoine and Rudik initial condition
Mi0(:,20)=[0.45 ; 0.55; 0; 0]*263*3.6644; %Golosov, calculated in file InitialValuesFAIR_GL_Golosov.m 20% of emissions flow in never decaying whereas 32% in decaying box with e-folding lifetime of 300 years
Mi0(:,21)=[0.373;0.314;0.311; 0]*263*3.6644; % Gerlagh&Liski, calculated in InitialValuesFAIR_GL_Golosov.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options= odeset('RelTol',1e-6);
for scenario=1:8
%JOOS GEOFFROY (and Lemoine&Rudik + Golosov decay) 
F_nonco2=interp1(Decennia2005,F_NONCO2(scenario+1,:),year);
F_nonco2_5y=interp1(Decennia2005,F_NONCO2(scenario+1,:),2015:5:2100);
M_Joos_LR_Golos_GL=zeros(21,86);
%DICE has a steadystate stock in the atmosphere of 588GtC and a 2015 stock of 851GtC with a difference of 263. 
for i=1:21
    [t,M]=ode45(@(t,M) odefcn_decay(t,M,Epath(scenario,:),year,J(i,:),1),[firstyear lastyear],Mi0(:,i),options);
    M=sum(transpose(M));
    Mtemporary=interp1(t,M(1,:),year);
    M_Joos_LR_Golos_GL(i,:)=transpose(Mtemporary);
end
M_Joosbestfit(scenario,:)=M_Joos_LR_Golos_GL(18,:);
M_LR(scenario,:)=M_Joos_LR_Golos_GL(19,:);
M_Golosov(scenario,:)=M_Joos_LR_Golos_GL(20,:);
M_GL(scenario,:)=M_Joos_LR_Golos_GL(21,:);
T_JoosGeof=zeros(16*16,86);
for k=1:16 %loop runs over carbon decay models
  for i=1:16 %loop runs over thermal inertia models
    Forcing=Climatesensitivity*G(i,3)*log((MequilibriumGtCO2+M_Joos_LR_Golos_GL(k,:))/MequilibriumGtCO2)/log(2)+F_nonco2;
    [t,T]=ode45(@(t,T) odefcn_inertia(t,T,Forcing,year,G(i,:)),[firstyear lastyear],[T0atm;T0ocean]);
    Ttemporary=interp1(t,T(:,1),year);
    T_JoosGeof(16*(k-1)+i,:)=transpose(Ttemporary);
  end 
end
T_sorted=sort(T_JoosGeof,1);
T_deciles=zeros(9,86);
for i=1:9
   T_deciles(i,:)=T_sorted(round(i*25.6),:);
end
%bestfit trajectory
Forcing=Climatesensitivity*G(17,3)*log((MequilibriumGtCO2+M_Joos_LR_Golos_GL(18,:))/MequilibriumGtCO2)/log(2)+F_nonco2;
[t,T]=ode45(@(t,T) odefcn_inertia(t,T,Forcing,year,G(17,:)),[firstyear lastyear],[T0atm;T0ocean],options);
Ttemporary=interp1(t,T(:,1),year);
T_JoosGeofbestfit(scenario,:)=transpose(Ttemporary);
%Lemoine&Rudik
Forcing=Climatesensitivity*G(19,3)*log((MequilibriumGtCO2+M_Joos_LR_Golos_GL(19,:))/MequilibriumGtCO2)/log(2)+F_nonco2;
[t,T]=ode45(@(t,T) odefcn_inertia(t,T,Forcing,year,G(19,:)),[firstyear lastyear],[T0atm;T0ocean],options);
Ttemporary=interp1(t,T(:,1),year);
T_LR(scenario,:)=transpose(Ttemporary);
%Golosov
T_Golosov(scenario,:)= 3.1 *log(1+M_Joos_LR_Golos_GL(20,:)/588/3.664)/log(2)+ 1/ G(18,3)*F_nonco2;
%Gerlagh & Liski
Forcing=Climatesensitivity*G(20,3)*log((MequilibriumGtCO2+M_Joos_LR_Golos_GL(21,:))/MequilibriumGtCO2)/log(2)+F_nonco2;
[t,T]=ode45(@(t,T) odefcn_inertia(t,T,Forcing,year,G(20,:)),[firstyear lastyear],[T0atm;T0ocean],options);
Ttemporary=interp1(t,T(:,1),year);
T_GL(scenario,:)=transpose(Ttemporary);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Gerlagh Liski Decay old
% a_GL=[0.163; 0.184; 0.449];eta_GL=[0; 0.074; 0.470];
% M_GL_boxes=zeros(3,10);
% M_GL_boxes(1:3,1)=[0.373;0.314;0.311]*263*3.6644; % obtained by running model on historical emissions 
% for i=1:9
%     M_GL_boxes(1:3,i+1)=a_GL*Epath(scenario,i*10-5)*10 + (ones(3,1)-eta_GL).*M_GL_boxes(1:3,i);
% end
% M_GL(scenario,:)=interp1([2015:10:2105],sum(M_GL_boxes),year,'spline');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DICE Decay
M_DICE_GtC=zeros(12,18);
M_DICE_GtC(1:3,1)=[851;460;1740];
b_Dice=[0.88 0.196 0; 0.12 0.797 0.001465; 0 0.007 0.998535];
for i=1:17
    M_DICE_GtC(1:3,i+1)=[5/3.6644*Epath(scenario,i*5-2) ;0;0] + b_Dice*M_DICE_GtC(1:3,i);
end
M_DICE(scenario,:)=interp1(firstyear:5:lastyear,M_DICE_GtC(1,:),year, 'spline')*3.6644;
%DICE2013 decay
%b12=0.088; b23=0.0025; MATeq=588; MUPeq=1350; MLOeq=10000; b_DICE_2013=[1-b12  b12*MATeq/MUPeq       0; b12  1-b12*MATeq/MUPeq-b23 b23*MUPeq/MLOeq;0   b23   1-b23*MUPeq/MLOeq];
M_DICE_GtC(10:12,1)=[851;1534;10010]; % adapted to have the same atmospheric concentration as DICE2016
b_Dice2013=[0.912,0.038328,0;0.088,0.959171,0.0003375;0,0.0025,0.9996625]; % from file GraphsDecayInertiaIRF2010.m
for i=1:17
    M_DICE_GtC(10:12,i+1)=[5/3.6644*Epath(scenario,i*5-2) ;0;0] + b_Dice2013*M_DICE_GtC(10:12,i);
end
M_DICE2013(scenario,:)=interp1(firstyear:5:lastyear,M_DICE_GtC(10,:),year, 'spline')*3.6644;

%DICE decay fitted to JOOS
b12=0.224; %Dice has 0.12
b23=0.092; %Dice has 0.007
MATeq=588;  
MUPeq=399; %Dice has 360
MLOeq=760; %MATeq/0.283-MUPeq; %to get 28.3% never-decaying proportion as in Dice
M_DICE_GtC(4:6,1)=[851;515;816];  %initial values recalculated on historic emissions, but MAT_2015 is 893 according to historic emissions, 3 stocks above preindustrial are reduced by the same %
b_Dice_Joosfit=[1-b12  b12*MATeq/MUPeq       0;...
                b12    1-b12*MATeq/MUPeq-b23 b23*MUPeq/MLOeq;...
                0      b23                   1-b23*MUPeq/MLOeq];
for i=1:17
    M_DICE_GtC(4:6,i+1)=[5/3.6644*Epath(scenario,i*5-2)*0.992 ;0;0] + b_Dice_Joosfit*M_DICE_GtC(4:6,i);
end
%DICE TEMP
T_DICE5y=zeros(2,18);
T_DICE5y(:,1)=[T0atm;0.0068];
Forcing=Climatesensitivity*1.1875*log(M_DICE_GtC(1,:)*3.6644/MequilibriumGtCO2)/log(2)+F_nonco2_5y;
for i=1:17
    T_DICE5y(1,i+1)=T_DICE5y(1,i)+0.1005*(Forcing(1,i+1)-1.1875*T_DICE5y(1,i)-0.088*(T_DICE5y(1,i)-T_DICE5y(2,i)));
    T_DICE5y(2,i+1)=T_DICE5y(2,i)+0.025*(T_DICE5y(1,i)-T_DICE5y(2,i));
end
T_DICE(scenario,:)=transpose(interp1(firstyear:5:lastyear,transpose(T_DICE5y(1,:)),year));
%DICE trajectory with Geoffroy parameters
T_DICE_Geoff=zeros(2,18);
T_DICE_Geoff(:,1)=[T0atm;T0ocean];
Forcing=Climatesensitivity*1.13*log(M_DICE_GtC(1,:)*3.664/MequilibriumGtCO2)/log(2)+F_nonco2_5y;
for i=1:17
    T_DICE_Geoff(1,i+1)=T_DICE_Geoff(1,i)+0.386*(Forcing(1,i+1)-1.13*T_DICE_Geoff(1,i)-0.73*(T_DICE_Geoff(1,i)-T_DICE_Geoff(2,i)));
    T_DICE_Geoff(2,i+1)=T_DICE_Geoff(2,i)+0.034*(T_DICE_Geoff(1,i)-T_DICE_Geoff(2,i));
end
T_DICE_Geoff=transpose(interp1(firstyear:5:lastyear,transpose(T_DICE_Geoff),year));
%DICE trajectory with Geoffroy parameters and decay fitted to Joos
T_DICE_JoosFitGeoff=[[T0atm;T0ocean] zeros(2,17)];
Forcing=Climatesensitivity*1.13*log(M_DICE_GtC(4,:)*3.664/MequilibriumGtCO2)/log(2)+F_nonco2_5y;
for i=1:17
    T_DICE_JoosFitGeoff(1,i+1)=T_DICE_JoosFitGeoff(1,i)+0.386*(Forcing(1,i+1)-1.13*T_DICE_JoosFitGeoff(1,i)-0.73*(T_DICE_JoosFitGeoff(1,i)-T_DICE_JoosFitGeoff(2,i)));
    T_DICE_JoosFitGeoff(2,i+1)=T_DICE_JoosFitGeoff(2,i)+0.034*(T_DICE_JoosFitGeoff(1,i)-T_DICE_JoosFitGeoff(2,i));
end
T_DICE_JoosFitGeoff=transpose(interp1(firstyear:5:lastyear,transpose(T_DICE_JoosFitGeoff),year));
%DICE fitted to FAIR
Ecum_GTC=[571 zeros(1,17)];
Alpha=zeros(1,18);
M_DICE_GtC(7:9,:)=[[851;515;816] zeros(3,17)];
T_DICE_FAIRFitGeoff=[[T0atm;T0ocean] zeros(2,17)];
for i=1:17
Ecum_GTC(i+1)=Ecum_GTC(i)+5*Epath(scenario,i*5-2)/3.664;
Alpha(i)=0.009653495*exp(0.08855*(34.4+0.019*(Ecum_GTC(i+1)-M_DICE_GtC(7,i)+MequilibriumGtCO2/3.664)+4.165*(T_DICE_FAIRFitGeoff(1,i))));
%Ecum_GTC(i) would be another option in line above
b12=0.224/Alpha(i); %Dice has 0.12
b23=0.092/Alpha(i);
b_Dice_FAIRfit=[1-b12 b12*MATeq/MUPeq       0;...
                b12    1-b12*MATeq/MUPeq-b23 b23*MUPeq/MLOeq;...
                0      b23                   1-b23*MUPeq/MLOeq];
M_DICE_GtC(7:9,i+1)=[5/3.6644*Epath(scenario,i*5-2)*0.992 ;0;0] + b_Dice_FAIRfit*M_DICE_GtC(7:9,i);
Forcing=Climatesensitivity*1.13*log(M_DICE_GtC(7,:)*3.6644/MequilibriumGtCO2)/log(2)+F_nonco2_5y;
T_DICE_FAIRFitGeoff(1,i+1)=T_DICE_FAIRFitGeoff(1,i)+0.386*(Forcing(1,i+1)-1.13*T_DICE_FAIRFitGeoff(1,i)-0.73*(T_DICE_FAIRFitGeoff(1,i)-T_DICE_FAIRFitGeoff(2,i)));
T_DICE_FAIRFitGeoff(2,i+1)=T_DICE_FAIRFitGeoff(2,i)+0.034*(T_DICE_FAIRFitGeoff(1,i)-T_DICE_FAIRFitGeoff(2,i));
end
T_DICE_FAIRFitGeoff=transpose(interp1(firstyear:5:lastyear,transpose(T_DICE_FAIRFitGeoff),year));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FAIR
%Alpha=transpose(interp1(Decennia, transpose(ALPHA) ,year));
[t,X]=ode45(@(t,X) odefcn_endogAlpha(t,X,year,Epath(scenario,:),F_nonco2,J(18,:),G(17,:),3.1),...
    [firstyear lastyear],[Mi0(:,1);T0atm;T0ocean;Ecum2015],options);
Ttemporary=interp1(t,X(:,5),year);
T_FAIR(scenario,:)=transpose(Ttemporary);
Mtemporary=sum(X(:,1:4),2);
M_FAIR(scenario,:)=transpose(interp1(t,Mtemporary,year));
X2100(scenario,:)=X(end,:);
end
%--------------------------------------------------------------------------------------------------------------------------------
%UNREPORTED PLOTS Emissions, Concentrations Temp per RCP scenario
subplot(3,1,1);
    plot(year,Epath(scenario,:),'b');
    title(txt_RCP(scenario+1));
    xlabel('years');
    ylabel('Emissions (GtCO2)');
    axis([2015 2100 0 inf]);
    axis 'auto y';
    pbaspect([3 1 1])
subplot(3,1,2);
    plot(year,transpose(M_Joos_LR_Golos_GL(1:16,:))*47/366+275.7,'color',[0 0.5 0.25],'linewidth',0.1);
    hold on
    p1=plot(year,transpose(M_Joos_LR_Golos_GL(1,:))*47/366+275.7,'color',[0 0.5 0.25],'linewidth',0.1);
    p2=plot(year,transpose(M_Joos_LR_Golos_GL(18,:))*47/366+275.7,'color',[0 0.2 0],'linewidth',1.5);
    p3=plot(year,interp1(firstyear:5:lastyear,(M_DICE_GtC(1,:))*47/100,year),'b','linewidth',1);
    p4=plot(year,interp1(firstyear:5:lastyear,(M_DICE_GtC(4,:))*47/100,year),'LineStyle','--','Color','b','LineWidth',0.5);
    p6=plot(year, M_FAIR(scenario,:)*47/366+275.7,'color', [1 0.7 0],'linewidth',1);
    p5=plot(year,interp1(firstyear:5:lastyear,(M_DICE_GtC(7,:))*47/100,year),'LineStyle','--','Color','b','LineWidth',0.5);
    p7=plot(year, interp1(Decennia2005,CONC(scenario+1,:),year),'color', [1 0 0],'linewidth',1);
    hold off
    lgd=legend([p1 p2 p3 p4 p5 p6 p7],'Joos','Best fit Joos et al.','DICE','DICE fitted to Joos','DICE fitted to FAIR',...
        'FAIR','MAGICC','Location','NorthWest');
    lgd.FontSize=4;
    xlabel('years')
    ylabel('Atmospheric CO2 (ppm)')
    axis([2015 2100 -inf inf]);
    pbaspect([3 1 1])
subplot(3,1,3);
    ar=area(transpose(year),transpose(T_deciles(9,:)), 'FaceColor' , [0.5 0.98 0.6]);
    hold on
    area(transpose(year),transpose(T_deciles(8,:)), 'FaceColor' , [0.5 0.96 0.6])
    area(transpose(year),transpose(T_deciles(7,:)), 'FaceColor' , [0.5 0.90 0.6])
    area(transpose(year),transpose(T_deciles(6,:)), 'FaceColor' , [0.5 0.80 0.6])
    area(transpose(year),transpose(T_deciles(4,:)), 'FaceColor' , [0.5 0.90 0.6])
    area(transpose(year),transpose(T_deciles(3,:)), 'FaceColor' , [0.5 0.96 0.6])
    area(transpose(year),transpose(T_deciles(2,:)), 'FaceColor' , [0.5 0.98 0.6])
    area(transpose(year),transpose(T_deciles(1,:)), 'FaceColor' , [1 1 1])
    p1=plot(year,T_JoosGeofbestfit(scenario,:),'color',[0 0.2 0],'linewidth',1.5);
    p2=plot(year, T_DICE(scenario,:),'b','linewidth',1);
    p3=plot(year, T_DICE_JoosFitGeoff(1,:),'LineStyle','--','Color','b','LineWidth',0.5);
    p4=plot(year, T_DICE_FAIRFitGeoff(1,:),'LineStyle','--','Color','b','LineWidth',0.5);
    p5=plot(year, T_FAIR(scenario,:),'color', [1 0.8 0],'linewidth',1);
    p6=plot(year, interp1(Decennia2005,TEMP(scenario+1,:),year),'color', [1 0 0],'linewidth',1);
    hold off
    lgd=legend([ar p1 p2 p3 p4 p5 p6],'Deciles Joos-Geoffroy combinations','Best fit Joos-Geoffroy ensemble',...
        'DICE','DICE with Geoffroy inertia and approximate Joos decay','DICE with Geoffroy inertia and approximate FAIR decay',...
        'FAIR with Geoffroy inertia','MAGICC','Location','SouthEast');
    lgd.FontSize=4;
    xlabel('years')
    ylabel('Global warming (°C)')
    axis([2015 2100 0 inf]);
    pbaspect([3 1 1])
print(txt_RCP{scenario+1,1},'-fillpage', '-dpdf' ); 
save('X2100_RCP.mat', 'X2100') ;


%Create variables for cumulative emissions
    M_FAIR_ppm=280 + M_FAIR*47/366.4;
    M_DICE_ppm=M_DICE*47/366.4;
    M_DICE2013_ppm=M_DICE2013*47/366.4;
    M_Joos_ppm=280 + M_Joosbestfit*47/366.4;
    M_FUND_ppm=[410.5372674	413.680392	416.8026076	419.9075855	422.9980366	426.0760556	429.1433321	432.2012804	435.2511203	438.2939279	441.3306676	444.3622135	447.3893638	450.4128503	453.4333458	456.45147	459.4677932	462.4828409	465.4970961	468.5110028	471.5249679	474.5393638	477.5545307	480.5707781	483.5883868	486.6076107	489.6286782	492.6517936	495.6771389	498.7048745	501.7351408	504.7680593	507.8037337	510.8422504	513.8836802	516.9280789	519.9754877	523.0259347	526.0794352	529.1359923	532.1955981	535.2582338	538.3238704	541.3924695	544.4639835	547.5383563	550.6155238	553.6954141	556.7779484	559.8630407	562.9505989	566.0405247	569.132714	572.2270575	575.3234407	578.4217444	581.5218448	584.623614	587.7269202	590.8316276	593.9375971	597.0446864	600.1527501	603.26164	606.3712051	609.4812922	612.5917458	615.7024081	618.8131196	621.9237191	625.0340435	628.1439285	631.2532086	634.3617168	637.4692852	640.5757453	643.6809274	646.7846613	649.8867764	652.9871015	656.0854652	659.1816958	662.2756215	665.3670706	668.4558712	671.541852];
    M_PAGE_ppm=[407.3396718	409.4699131	411.5852758	413.6860377	415.7724723	417.844849	419.9034333	421.9484865	423.9802662	425.9990258	428.0050151	429.9984799	431.9796624	433.9488011	435.9061306	437.851882	439.786283	441.7095575	443.6219049	445.5235212	447.4145989	449.2953274	451.1658927	453.0264776	454.877262	456.7184224	458.5501325	460.3725627	462.1858809	463.9902517	465.7858371	467.5727962	469.3512853	471.1214583	472.8834659	474.6374567	476.3835765	478.1219684	479.8527734	481.5761297	483.2921732	485.0010376	486.702854	488.3977513	490.0858562	491.7672933	493.4421847	495.1106506	496.7728091	498.4287761	500.0786656	501.7225894	503.3606576	504.9929782	506.6196573	508.2407992	509.8565063	511.4668793	513.072017	514.6720165	516.2669732	517.8569809	519.4421317	521.0225159	522.5982225	524.1693388	525.7359504	527.2981417	528.8559953	530.4095927	531.9590136	533.5043365	535.0456385	536.5829953	538.1164813	539.6461695	541.1721317	542.6944385	544.2131592	545.7283617	547.240113	548.7484788	550.2535235	551.7553107	553.2539025	554.7493602];
    M_LR_ppm=280 + M_LR*47/366.4;
    M_Golosov_ppm=280 + M_Golosov*47/366.4;
    M_GL_ppm=280+M_GL*47/366.4;
    %cumulative absorbed emissions
    M_FAIR_uptake_cum=Ecum_yearly - M_FAIR;
    M_DICE_uptake_cum=Ecum_yearly - M_DICE + MequilibriumGtCO2 ;
    M_DICE2013_uptake_cum=Ecum_yearly - M_DICE2013 + MequilibriumGtCO2 ;
    M_Joos_uptake_cum=Ecum_yearly - M_Joosbestfit;
    M_FUND_uptake_cum=[1063.065242	1077.670165	1092.438388	1107.341237	1122.35754	1137.470937	1152.668234	1167.938385	1183.271861	1198.66026	1214.09605	1229.572402	1245.083084	1260.622381	1276.185038	1291.766214	1307.361456	1322.96666	1338.578053	1354.192168	1369.805827	1385.416121	1401.020394	1416.616228	1432.20143	1447.774018	1463.332208	1478.874402	1494.399182	1509.905293	1525.39164	1540.857272	1556.301382	1571.723293	1587.122452	1602.498425	1617.850888	1633.179623	1648.484511	1663.765526	1679.022731	1694.256273	1709.466377	1724.653345	1739.817547	1754.959421	1770.079469	1785.178252	1800.256386	1815.314541	1830.353438	1845.373844	1860.376572	1875.362475	1890.332449	1905.287423	1920.228365	1935.156274	1950.072179	1964.977141	1979.872245	1994.758604	2009.637352	2024.509649	2039.376671	2054.239617	2069.099701	2083.958155	2098.816224	2113.675169	2128.536261	2143.400786	2158.270035	2173.145312	2188.027929	2202.919202	2217.820456	2232.73302	2247.658228	2262.597415	2277.551921	2292.523086	2307.512253	2322.520762	2337.549955	2352.601171];
    M_PAGE_uptake_y=[1111.464495	1133.976069	1156.603844	1179.345652	1202.199356	1225.162851	1248.234066	1271.410958	1294.691517	1318.073763	1341.555745	1365.135542	1388.811265	1412.58105	1436.443065	1460.395504	1484.43659	1508.564574	1532.777899	1557.075034	1581.454475	1605.914744	1630.454387	1655.071976	1679.766108	1704.535403	1729.378506	1754.294084	1779.280827	1804.337449	1829.462685	1854.655293	1879.91405	1905.237758	1930.625236	1956.075326	1981.586889	2007.158805	2032.789977	2058.479322	2084.22578	2110.028307	2135.88588	2161.79749	2187.762148	2213.778883	2239.846739	2265.964778	2292.132078	2318.347734	2344.610855	2370.920568	2397.276013	2423.676347	2450.120742	2476.608382	2503.138467	2529.710212	2556.322845	2582.975607	2609.667753	2636.39855	2663.167281	2689.973238	2716.815728	2743.694068	2770.607589	2797.555633	2824.537554	2851.552717	2878.600497	2905.680284	2932.791473	2959.933475	2987.105707	3014.3076	3041.538593	3068.798135	3096.085685	3123.400711	3150.742691	3178.111112	3205.505471	3232.925272	3260.370027	3287.839261];
    M_LR_uptake_cum=Ecum_yearly - M_LR;
    M_Golosov_uptake_cum=Ecum_yearly - M_Golosov;
    M_GL_uptake_cum=Ecum_yearly - M_GL;
    %yearly uptake
    M_FAIR_uptake_y= Epath(:,2:86)-(M_FAIR(:,2:86)-M_FAIR(:,1:85));
    M_DICE_uptake_y= Epath(:,2:86)-(M_DICE(:,2:86)-M_DICE(:,1:85));
    M_Joos_uptake_y=Epath(:,2:86) -(M_Joosbestfit(:,2:86)-M_Joosbestfit(:,1:85));
    M_FUND_uptake_y=[13.88371409	14.6049232	14.76822247	14.90284875	15.01630314	15.11339775	15.19729696	15.27015013	15.33347643	15.38839927	15.43578931	15.4763522	15.51068248	15.53929706	15.56265622	15.58117694	15.59524163	15.60520394	15.61139279	15.61411533	15.61365909	15.6102938	15.6042728	15.59583426	15.58520226	15.57258778	15.55818951	15.54219469	15.52477981	15.50611131	15.48634614	15.46563243	15.44410994	15.42191063	15.39915908	15.37597293	15.35246332	15.32873522	15.30488784	15.28101489	15.25720496	15.23354175	15.21010437	15.18696759	15.16420207	15.14187457	15.12004818	15.09878249	15.07813378	15.0581552	15.0388969	15.02040622	15.00272777	14.98590364	14.96997344	14.95497448	14.94094185	14.92790852	14.91590543	14.9049616	14.89510423	14.88635871	14.87874877	14.87229653	14.86702252	14.86294583	14.86008409	14.85845357	14.85806923	14.85894474	14.86109258	14.86452404	14.86924928	14.87527739	14.8826164	14.89127333	14.90125421	14.91256415	14.92520735	14.93918713	14.95450594	14.97116545	14.98916652	15.00850923	15.02919294	15.05121629];
    M_PAGE_uptake_y=[22.09653839	22.51157353	22.62777512	22.7418078	22.85370401	22.96349579	23.07121478	23.17689222	23.28055894	23.38224541	23.48198174	23.57979764	23.6757225	23.76978531	23.86201474	23.9524391	24.04108636	24.12798413	24.2133248	24.29713503	24.37944103	24.46026858	24.53964304	24.61758934	24.69413203	24.76929523	24.84310267	24.91557771	24.98674331	25.05662205	25.12523617	25.1926075	25.25875757	25.3237075	25.38747809	25.45008982	25.51156278	25.57191678	25.63117127	25.68934539	25.74645796	25.80252749	25.85757219	25.91160994	25.96465835	26.01673471	26.06785604	26.11803906	26.16730022	26.21565568	26.26312133	26.30971278	26.3554454	26.40033428	26.44439424	26.48763986	26.53008547	26.57174516	26.61263275	26.65276185	26.69214582	26.73079778	26.76873063	26.80595705	26.84248949	26.87834018	26.91352115	26.94804419	26.98192092	27.01516271	27.04778077	27.07978609	27.11118946	27.1420015	27.17223261	27.20189304	27.23099282	27.25954184	27.28754977	27.31502615	27.34198031	27.36842143	27.39435854	27.41980048	27.44475595	27.46923347];
    M_LR_uptake_y= Epath(:,2:86)-(M_LR(:,2:86)-M_LR(:,1:85));
    M_LR_uptake_y(8,:)=interp1(M_LR_ppm(8,1:14:85),M_LR_uptake_y(8,1:14:85),M_LR_ppm(8,1:85),'spline');
    M_Golosov_uptake_y= Epath(:,2:86)-(M_Golosov(:,2:86)-M_Golosov(:,1:85));
    M_GL_uptake_y= Epath(:,2:86)-(M_GL(:,2:86)-M_GL(:,1:85));

%FIG 2: Plot yearly uptake by carbon sinks against concentration
%constant 2015 emissions
    subplot(1,1,1);
    p1=plot(M_FAIR_ppm(8,1:85),M_FAIR_uptake_y(8,:),'r','linewidth',3);
        hold on
    scatter(M_FAIR_ppm(8,1:5:85),M_FAIR_uptake_y(8,1:5:85),'*','r','linewidth',2);
    plot(M_DICE_ppm(8,1:85),M_DICE_uptake_y(8,:),'k','linewidth',1.5) ;
        p3=scatter(M_DICE_ppm(8,1:5:85),M_DICE_uptake_y(8,1:5:85),'s','MarkerEdgeColor',[0.25 0.6 1],'linewidth',1.5) ;
    plot(M_DICE2013_ppm(8,1:85),M_DICE_uptake_y(8,:),'k','linewidth',1.5) ;
        p4=scatter(M_DICE2013_ppm(8,1:5:85),M_DICE_uptake_y(8,1:5:85),'p','MarkerEdgeColor','c','linewidth',1.5) ;
    plot(M_FUND_ppm,M_FUND_uptake_y,'k','linewidth',1.5) ;
        p5=scatter(M_FUND_ppm(1,1:5:85),M_FUND_uptake_y(1,1:5:85),'v','MarkerEdgeColor','m','linewidth',1.5) ;
    plot(M_PAGE_ppm,M_PAGE_uptake_y,'k','linewidth',1.5) ;
        p6=scatter(M_PAGE_ppm(1,1:5:85),M_PAGE_uptake_y(1,1:5:85),'d','MarkerEdgeColor','r','linewidth',1.5) ;
    plot(M_Golosov_ppm(8,1:85),M_Golosov_uptake_y(8,1:85),'k','linewidth',1.5)
        p7=scatter(M_Golosov_ppm(8,1:5:85),M_Golosov_uptake_y(8,1:5:85),'+','MarkerEdgeColor',[0.9 0.55 0],'linewidth',1.5) ;
    plot(M_GL_ppm(8,1:5:85),M_GL_uptake_y(8,1:5:85),'k','linewidth',1.5) ;
        p8=scatter(M_GL_ppm(8,1:5:85),M_GL_uptake_y(8,1:5:85),'x','MarkerEdgeColor',[1 0.9 0.1],'linewidth',1.5) ;

    plot(M_LR_ppm(8,1:5:85), M_LR_uptake_y(8,1:5:85),'k','linewidth',1.5) ;
        p9=scatter(M_LR_ppm(8,1:5:85), M_LR_uptake_y(8,1:5:85),'^','MarkerEdgeColor','b','linewidth',1.5) ;
    %plot(M_LR_ppm(8,1:14:end), interp1(M_LR_ppm(8,1:85),M_LR_uptake_y(8,:),M_LR_ppm(8,1:14:end)),'-.r','linewidth',1.5) ;
    hold off
    lgd=legend([p1 p3 p4 p5 p6 p7 p8 p9],'FAIR','DICE16','DICE13','FUND','PAGE','GHKT14','GL18','LR17',... 
        'FontSize',11,'NumColumns',2,'Location','Southeast');
    ylabel('Yearly Emissions absorbed by sinks (GtCO2)')
    xlabel('CO2 in atmosphere (ppm)')
    axis([400 700 5 35]);
    set(gcf,'PaperOrientation','landscape');
    print('-fillpage','Fig2', '-dpdf')
  
%FIG Appendix E: Warming against cumulative emissions
%IPCC vs DICE FUND PAGE FAIR
Labels={'RCP 2.6','RCP 2.6','RCP 4.5','','RCP 6.0','RCP 6.0','' };
subplot(1,1,1);
    for i=[1 3 5 7] ; p1=plot(IPCC_cum(i+1,:),IPCC_cum(i,:),'b');t=text(max(IPCC_cum(i+1,:)),max(IPCC_cum(i,:)),Labels(i));t.FontSize=8;hold on; end
    for i=[2 3 6 7] ; p2=plot(Ecum_yearly(i,1:81),T_DICE(i,1:81),':r','LineWidth',2) ;t=text(max(Ecum_yearly(i,1:81)),max(T_DICE(i,1:81)),Labels(i));t.FontSize=8; end
    t=text(7400,4.3,'RCP 8.5');t.FontSize=8;
    t=text(7400,3.9,'RCP 8.5');t.FontSize=8;
    hold off
    legend([p1 p2],'IPCC AR5 model mean','DICE','Location','Southeast');
    xlabel('Cumulative Emissions GtCO2')
    ylabel('Global warming (°C)')
    axis([2000,8000,0.7,4.5])
    set(gcf,'PaperOrientation','landscape');   
    pbaspect([1 0.7 1])
    print('-fillpage','CumEmissions_vs_Temp_DICE', '-dpdf');   
subplot(1,1,1);
    for i=[1 3 5 7] ; p1=plot(IPCC_cum(i+1,:),IPCC_cum(i,:),'b');t=text(max(IPCC_cum(i+1,:)),max(IPCC_cum(i,:)),Labels(i));t.FontSize=8;hold on; end
    for i=[1 3 5 7] ; p2=plot(FUND_cum(i+1,1:81),FUND_cum(i,1:81),':r','LineWidth',2) ;t=text(max(FUND_cum(i+1,1:81)),max(FUND_cum(i,1:81)),Labels(i));t.FontSize=8; end
    t=text(7400,4.3,'RCP 8.5');t.FontSize=8;
    t=text(7400,3.9,'RCP 8.5');t.FontSize=8;
    hold off
    legend([p1 p2],'IPCC AR5 model mean','FUND','Location','Southeast');
    xlabel('Cumulative Emissions GtCO2')
    ylabel('Global warming (°C)')
    axis([2000,8000,0.7,4.5])
    set(gcf,'PaperOrientation','landscape');
    pbaspect([1 0.7 1])
    print('-fillpage','CumEmissions_vs_Temp_FUND', '-dpdf');   
subplot(1,1,1);
    for i=[1 3 5 7] ; p1=plot(IPCC_cum(i+1,:),IPCC_cum(i,:),'b');t=text(max(IPCC_cum(i+1,:)),max(IPCC_cum(i,:)),Labels(i));t.FontSize=8; hold on; end
    for i=[1 3 5 7] ; p2=plot(PAGE_cum(i+1,1:81),PAGE_cum(i,1:81),':r','LineWidth',2) ;t=text(max(PAGE_cum(i+1,1:81)),max(PAGE_cum(i,1:81)),Labels(i));t.FontSize=8; end
    t=text(7400,2.4,'RCP 8,5');t.FontSize=8;
    t=text(7400,3.9,'RCP 8,5');t.FontSize=8;
    hold off
    legend([p1 p2],'IPCC AR5 model mean','PAGE','Location','Southeast');
    xlabel('Cumulative Emissions GtCO2')
    ylabel('Global warming (°C)')
    axis([2000,8000,0.7,4.5])
    set(gcf,'PaperOrientation','landscape');
    pbaspect([1 0.7 1])
    print('-fillpage','CumEmissions_vs_Temp_PAGE', '-dpdf');   
subplot(1,1,1);
    for i=[1 3 5 7] ; p1=plot(IPCC_cum(i+1,:),IPCC_cum(i,:),'b');t=text(max(IPCC_cum(i+1,:)),max(IPCC_cum(i,:)),Labels(i));t.FontSize=8; hold on; end
    for i=[2 3 6 7] ; p2=plot(Ecum_yearly(i,1:81),T_LR(i,1:81),':r','LineWidth',2) ;t=text(max(Ecum_yearly(i,1:81)),max(T_LR(i,1:81)),Labels(i));t.FontSize=8; end
    t=text(7400,2.5,'RCP 8.5');t.FontSize=8;
    t=text(7400,3.9,'RCP 8.5');t.FontSize=8;
    hold off
    legend([p1 p2],'IPCC AR5 model mean','Lemoine&Rudik','Location','Southeast');
    xlabel('Cumulative Emissions GtCO2')
    ylabel('Global warming (°C)')
    axis([2000,8000,0.7,4.5])
    set(gcf,'PaperOrientation','landscape');
    pbaspect([1 0.7 1])
    print('-fillpage','CumEmissions_vs_Temp_LR', '-dpdf');   
subplot(1,1,1);
    for i=[1 3 5 7] ; p1=plot(IPCC_cum(i+1,:),IPCC_cum(i,:),'b');t=text(max(IPCC_cum(i+1,:)),max(IPCC_cum(i,:)),Labels(i));t.FontSize=8; hold on; end
    for i=[2 3 6 7] ; p2=plot(Ecum_yearly(i,1:81),T_FAIR(i,1:81),':r','LineWidth',2) ;t=text(max(Ecum_yearly(i,1:81)),max(T_FAIR(i,1:81)),Labels(i));t.FontSize=8; end
    t=text(7400,4.3,'RCP 8.5');t.FontSize=8;
    t=text(7400,3.9,'RCP 8.5');t.FontSize=8;
    hold off
    legend([p1 p2],'IPCC AR5 model mean','FAIR','Location','Southeast');
    xlabel('Cumulative Emissions GtCO2')
    ylabel('Global warming (°C)')
    axis([2000,8000,0.7,4.5])
    set(gcf,'PaperOrientation','landscape');
    pbaspect([1 0.7 1])
    print('-fillpage','CumEmissions_vs_Temp_FAIR', '-dpdf');   
subplot(1,1,1);
    for i=[1 3 5 7] ; p1=plot(IPCC_cum(i+1,:),IPCC_cum(i,:),'b');t=text(max(IPCC_cum(i+1,:)),max(IPCC_cum(i,:)),Labels(i));t.FontSize=8; hold on;end
    t=text(7400,3,Labels(7));t.FontSize=8;
    for i=[2 3 6 7] ; p2=plot(Ecum_yearly(i,1:81),T_Golosov(i,1:81),':r','LineWidth',2) ;t=text(max(Ecum_yearly(i,1:81)),max(T_Golosov(i,1:81)),Labels(i));t.FontSize=8; end
    t=text(7400,5.8,'RCP 8.5');t.FontSize=8;
    t=text(7400,3.9,'RCP 8.5');t.FontSize=8;
    hold off
    legend([p1 p2],'IPCC AR5 model mean','GHKT14','Location','Southeast');
    xlabel('Cumulative Emissions GtCO2')
    ylabel('Global warming (°C)')
    axis([2000,8000,0.7,6])
    set(gcf,'PaperOrientation','landscape');
    pbaspect([1 0.7 1])
    print('-fillpage','CumEmissions_vs_Temp_Golosov', '-dpdf');   
subplot(1,1,1);
    for i=[1 3 5 7] ; p1=plot(IPCC_cum(i+1,:),IPCC_cum(i,:),'b');t=text(max(IPCC_cum(i+1,:)),max(IPCC_cum(i,:)),Labels(i));t.FontSize=8; hold on;end
    for i=[2 3 6 7] ; p2=plot(Ecum_yearly(i,1:81),T_GL(i,1:81),':r','LineWidth',2) ;t=text(max(Ecum_yearly(i,1:81)),max(T_GL(i,1:81)),Labels(i));t.FontSize=8; end
    t=text(7600,3.2,'RCP 8.5');t.FontSize=8;
    t=text(7600,4.3,'RCP 8.5');t.FontSize=8;
    hold off
    legend([p1 p2],'IPCC AR5 model mean','GL18','Location','Southeast');
    xlabel('Cumulative Emissions GtCO2')
    ylabel('Global warming (°C)')
    axis([2000,8000,0.7,4.5])
    set(gcf,'PaperOrientation','landscape');
    pbaspect([1 0.7 1])
    print('-fillpage','CumEmissions_vs_Temp_GerlaghLiski', '-dpdf');     
