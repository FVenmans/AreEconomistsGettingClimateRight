%THIS FILE MAKES FIG 1, 3 AND 4 IN THE PAPER
%IMPACT RESPONSE FUNCTION FOR A PULSE OF EMISSIONS OF 100GtC IN 2010 
%To make pulse of 1 or 1000GtC: see file Fig1_PulseSize_Moment
close all; clear;
cd 'C:\Users\530443\Google Drive\Research Projects\Critique DICE\GraphsDecayInertia';% change here for the path
%M is always GtCO2
%READ XLS FILES
    [J,txt_J,Table_J] = xlsread('decay and inertia parameters.xlsx',1,'A1:H19','basic'); %import parameters from Joos et al.
    [G,txt_G,Table_G] = xlsread('decay and inertia parameters.xlsx',3,'A2:G20');%import parameters for Geoffroy et al.
    G(19,:)=G(18,:);G(18,:)=G(17,:);G(17,:)=G(19,:); %set best fit on line 18, as in J
    load('RCP.mat')%load background data of RCP scenarios
    T_PAGE=xlsread('Page Cycle SD Sept 2019.xlsx',2,'C5:GU5','basic');
    T_FUND_Beta2600=xlsread('FUND 3.11.xlsx',3,'C17:GU17','basic');
    M_FUND_diff_ppm=xlsread('FUND and PAGE output for figures 3 and 4.xlsx',1,'D6:GV6','basic');
    T_cf_FUND_diff=xlsread('FUND and PAGE output for figures 3 and 4.xlsx',1,'C12:GU12','basic');
    M_PAGE_diff_ppm=xlsread('FUND and PAGE output for figures 3 and 4.xlsx',2,'D7:GV7','basic')*46.9 ; % the excell line has weird unit, 29000 at the start
    T_cf_PAGE_diff=xlsread('FUND and PAGE output for figures 3 and 4.xlsx',2,'C12:GU12','basic');
%PARAMETERS
    Pulse=100*3.664; %DICE has 399.4 ppm and 851GtC in 2015=> Pulse of 100GtC=366.3GtCO2=46,9ppm    
    Mequilibrium=588*3.664; %Dice equilibrium concentration = 588GtC*0.469ppm/GtC
    FCO22x=3.6813; %forcing for a doubling of CO2 concentration from DICE
    Climatesensitivity=3.1; %value in DICE
    year=0:1:200 ;
%INITIAL CONDITIONS
    T0atm=0.85; %Atmosphere warming in 2015 in DICE is 0.85
    T0ocean=0.22; 
    T0ocean_DICE=0.0068; %Deep ocean warming in 2015 in DICE is 0.0068
    M2010=389/46.9*366.4; %background concentration in GtCO2 before pulse (from Joos 389ppm)
    Ecum2010=531*3.664; %value from Fair excel file. Richard Millar in Nature has 545GtC in 2015
    Mi0_DiceGtC=[851;460;1740];
    Mi0_Dice2013GtC=[830;1527;10010]; 
    Mi0=[0.473184713605214;0.359104481245797;0.143649591297395;0.0240612138515941]*(M2010-Mequilibrium);% FAIR calculated from historic emissions in InitialValuesFAIR_GL_Golosov
    %Mi0 =[0.5294;0.3433;0.1109;0.01622]*(M2010-Mequilibrium); % proportions received from Richard Millar, 
%EMISSIONS WHICH KEEP CONCENTRATION CONSTANT for median parameters  
    [t,M]=ode45(@(t,M) odefcn_decay_constM(t,M,J(18,:)),[0 200],Mi0);
    Mtemporary=interp1(t,M,year)';
    Epath=J(18,5:7)*Mtemporary(2:4,:);
    %Epath=zeros(1,201); %better for pulse in 2100 on RCP4.5
% %INITIAL CONDITIONS FOR RCP7(SSP3baseline) IN 2100 (line 4 in X2100) 
%     Mi0=[1643.6;1450.7;983.7;270.3]; 
%     M2010=Mequilibrium+sum(Mi0);  
%     T0atm=4.06;
%     T0ocean=1.36; 
%     T0ocean_DICE=1.36; 
%     Ecum2010=7565; 
%     Mi0_DiceGtC=(Ecum2010/3.664-531)*[0.63;0.2;0.17]+[588;360;1720];  
%     Mi0_Dice2013GtC=Mi0_DiceGtC;
% %INITIAL CONDITIONS FOR RCP4.5SSP2 IN 2100 (line 3 in X2100) 
%     Mi0=[1096.88;880.1;356.14;29.85]; 
%     M2010=Mequilibrium+sum(Mi0); 
%     T0atm=2.67;
%     T0ocean=1.12; 
%     T0ocean_DICE=1.12; 
%     Ecum2010=5045; 
%     Mi0_DiceGtC=(Ecum2010/3.664-531)*[0.63;0.2;0.17]+[588;360;1720];  %distribution over 3 boxes fits atmospheric carbon, shallow and deep ocean comparable to 2010
%     Mi0_Dice2013GtC=Mi0_DiceGtC;
%-------------------------------------------------------------------------------------------------------
%CARBON ABSORPTION (decay) (shown in Fig 3)
%JOOS et al. + Lemoine&Rudik + Golosov
M_JoosLRGolosGL_pulse=zeros(21,201);
M_JoosLRGolosGL_bg=zeros(21,201);
J(19,:)=[0 1 0 0 0.0138 0 0]; %Lemoine and Rudik
J(20,:)=[0.2 0.314 0.486 0	0.00228 10	0]; %Golosov
J(21,:)=[0.163 0.184 0.449 0 -0.1*log(1-0.074) -0.1*log(1-0.47) 0]; %Gerlagh&Liski delta=-1/10*ln(1-eta)
options= odeset('RelTol',1e-7);
for i=1:21
    %decay with pulse
    Mi0_pulse=Mi0+[J(i,1);J(i,2);J(i,3);J(i,4)]*Pulse; 
    [t,M]=ode45(@(t,M) odefcn_decay(t,M,Epath,year,J(i,:),1),[0 200],Mi0_pulse,options);
    M_JoosLRGolosGL_pulse(i,:)=interp1(t,sum(M,2),year);
    %background decay    
    [t,M]=ode45(@(t,M) odefcn_decay(t,M,Epath,year,J(i,:),1),[0 200],Mi0,options);
    M_JoosLRGolosGL_bg(i,:)=interp1(t,sum(M,2),year);
end
M_JoosLRGolosGL_diff=M_JoosLRGolosGL_pulse - M_JoosLRGolosGL_bg;
M_JoosLRGolosGL_diff(20,1)=M_JoosLRGolosGL_diff(20,2); %this corrects for immediate decay at time zero in Golosov et al.
%Initial values are incorrect for LR and Golosov, but difference is unaffected by initial stocks. The following gives the same result: 
Mi0_pulse=[J(20,1);J(20,2);J(20,3);J(20,4)]*Pulse; 
[t,M]=ode45(@(t,M) odefcn_decay(t,M,Epath,year,J(20,:),1),[0 200],Mi0_pulse,options); 
M_Golosov_diff(1,:)=interp1(t,sum(M,2),year);

%FAIR Decay and Temperature
%2010 pulse
T_FAIR_pulse=zeros(7,201);
T_FAIR_bg=zeros(7,201);
F_nonco2=zeros(1,201);
Xi0_pulse=[Mi0+[J(18,1);J(18,2);J(18,3);J(18,4)]*Pulse;T0atm;T0ocean;Ecum2010+Pulse];
Xi0_bg=[Mi0;T0atm;T0ocean;Ecum2010];
[t,X]=ode45(@(t,X) odefcn_endogAlpha(t,X,year,Epath,F_nonco2,J(18,:),G(18,:),3.1),[0 200],Xi0_pulse,options);
M_FAIR_pulse(1,:)=interp1(t,sum(X(:,1:4),2),year);
T_FAIR_pulse(1,:)=interp1(t,X(:,5),year);
[t,X]=ode45(@(t,X) odefcn_endogAlpha(t,X,year,Epath,F_nonco2,J(18,:),G(18,:),3.1),[0 200],Xi0_bg,options);
M_FAIR_bg(1,:)=interp1(t,sum(X(:,1:4),2),year);
T_FAIR_bg(1,:)=interp1(t,X(:,5),year);
%2100 pulse
load('X2100_RCP.mat')
X2100_pulse=X2100+ones(8,1)*[J(18,1) J(18,2) J(18,3) J(18,4) 0 0 1]*Pulse;
for k=2:8 %last scenario probably has a problematic alpha    
    [t,X]=ode45(@(t,X) odefcn_endogAlpha(t,X,year,Epath,F_nonco2,J(18,:),G(18,:),3.1),[0 200],X2100_pulse(k-1,:),options);
    M_FAIR_pulse(k,:)=interp1(t,sum(X(:,1:4),2),year);
    T_FAIR_pulse(k,:)=interp1(t,X(:,5),year);
    %background scenario
    [t,X]=ode45(@(t,X) odefcn_endogAlpha(t,X,year,Epath,F_nonco2,J(18,:),G(18,:),3.1),[0 200],X2100(k-1,:),options);
    M_FAIR_bg(k,:)=interp1(t,sum(X(:,1:4),2),year);
    T_FAIR_bg(k,:)=interp1(t,X(:,5),year);
end
M_FAIR_diff=M_FAIR_pulse-M_FAIR_bg    ;
T_FAIR_diff=T_FAIR_pulse-T_FAIR_bg    ;

%DICE Decay
M_DICE_GtC_pulse=ones(3,41);
M_DICE_GtC_pulse(:,1)=Mi0_DiceGtC +[Pulse/3.664;0;0];% [851+Pulse/3.664;460;1740];  %values 2015 in DICE.Change to 2010 values would not change the impact response function
b_DICE=[0.88 0.196 0; 0.12 0.797 0.001465; 0 0.007 0.998535];
for i=1:40
    M_DICE_GtC_pulse(:,i+1)=[5/3.664*Epath(5*i);0;0]+ b_DICE*M_DICE_GtC_pulse(:,i);
end
M_DICE_pulse(1,:)=interp1(0:5:200,M_DICE_GtC_pulse(1,:),year)*3.664; 
%background decay
M_DICE_GtC_bg=ones(3,41);
M_DICE_GtC_bg(:,1)=Mi0_DiceGtC;  
for i=1:40
    M_DICE_GtC_bg(:,i+1)= [5/3.664*Epath(5*i);0;0]+ b_DICE*M_DICE_GtC_bg(:,i);
end
M_DICE_bg(1,:)=interp1(0:5:200,M_DICE_GtC_bg(1,:),year)*3.664; 
M_DICE_diff=M_DICE_pulse - M_DICE_bg;

%DICE2013 Decay
b12=0.088; b23=0.0025; MATeq=588; MUPeq=1350; MLOeq=10000; 
b_DICE_2013=[1-b12  b12*MATeq/MUPeq       0; b12  1-b12*MATeq/MUPeq-b23 b23*MUPeq/MLOeq;0   b23   1-b23*MUPeq/MLOeq];
M_DICE2013_GtC_pulse=ones(3,41);
M_DICE2013_GtC_pulse(:,1)=Mi0_Dice2013GtC+[Pulse/3.664;0;0];   %[830+Pulse/3.664;1527;10010];  %values 2015 in DICE
for i=1:40
    M_DICE2013_GtC_pulse(:,i+1)=[5/3.664*Epath(5*i);0;0]+ b_DICE_2013*M_DICE2013_GtC_pulse(:,i);
end
M_DICE2013_pulse(1,:)=interp1(0:5:200,M_DICE2013_GtC_pulse(1,:),year)*3.664; 
%background decay
M_DICE2013_GtC_bg=ones(3,41);
M_DICE2013_GtC_bg(:,1)=Mi0_Dice2013GtC; %[830;1527;10010];  
for i=1:40
    M_DICE2013_GtC_bg(:,i+1)=[5/3.664*Epath(5*i);0;0]+ b_DICE_2013*M_DICE2013_GtC_bg(:,i);
end 
M_DICE2013_bg(1,:)=interp1(0:5:200,M_DICE2013_GtC_bg(1,:),year)*3.664; 
M_DICE2013_diff=(M_DICE2013_pulse - M_DICE2013_bg);
%--------------------------------------------------------------------------------------------------------------
%THERMAL INERTIA FOR CONSTANT FORCING (shown further in figure 4)
%GEOFFROY et al. /Lemoine&Rudik /GL18
G(19,:)=[1 1 0.0091 0 0 0]; ... %Lemoine and Rudik parameters
%             1 1 0.0182 0 0 0; ... %Lemoine and Rudik high and low inertia
%             1 1 0.0046 0 0 0];
G(20,:)=[0.03 0.05 1.13 0 0 0]; % lambda is irrelevant, but should not be zero.
G(21,:)=[1 1 -0.1*log(1-0.183) 0 0 0];%Gerlagh and Liski epsilon=-1/10ln(1-0.183)
T_cf_JoosGeof_LR_GL_pulse=zeros(21,201);
for i=1:21
Forcing=Climatesensitivity*G(i,3)*log((M2010+Pulse)/Mequilibrium)/log(2)*ones(1,201);
[t,T]=ode45(@(t,T) odefcn_inertia(t,T,Forcing,year,G(i,:)),[0 200],[T0atm;T0ocean],options);
T_cf_JoosGeof_LR_GL_pulse(i,:)=transpose(interp1(t,T(:,1),year));
end 
%background temperature
T_cf_JoosGeof_LR_GL_bg=zeros(21,201);
for i=1:21
Forcing=Climatesensitivity*G(i,3)*log((M2010)/Mequilibrium)/log(2)*ones(1,201);
[t,T]=ode45(@(t,T) odefcn_inertia(t,T,Forcing,year,G(i,:)),[0 200],[T0atm;T0ocean],options);
T_cf_JoosGeof_LR_GL_bg(i,:)=transpose(interp1(t,T(:,1),year));
end 
T_cf_JoosGeof_LR_GL_diff=T_cf_JoosGeof_LR_GL_pulse - T_cf_JoosGeof_LR_GL_bg;

%DICE 
T_cf_DICE5y_pulse=[[T0atm;T0ocean_DICE] zeros(2,40)];
Fpulse=FCO22x*log((M2010+Pulse)/Mequilibrium)/log(2);
for i=1:40
    T_cf_DICE5y_pulse(1,i+1)=T_cf_DICE5y_pulse(1,i)+0.1005*(Fpulse-1.1875*T_cf_DICE5y_pulse(1,i)-0.088*(T_cf_DICE5y_pulse(1,i)-T_cf_DICE5y_pulse(2,i)));
    T_cf_DICE5y_pulse(2,i+1)=T_cf_DICE5y_pulse(2,i)+0.025*(T_cf_DICE5y_pulse(1,i)-T_cf_DICE5y_pulse(2,i));
end
%DICE Background temp
T_cf_DICE5y_bg=[[T0atm;T0ocean_DICE] zeros(2,40)];
Fpulse=FCO22x*log(M2010/Mequilibrium)/log(2);
for i=1:40
    T_cf_DICE5y_bg(1,i+1)=T_cf_DICE5y_bg(1,i)+0.1005*(Fpulse-1.1875*T_cf_DICE5y_bg(1,i)-0.088*(T_cf_DICE5y_bg(1,i)-T_cf_DICE5y_bg(2,i)));
    T_cf_DICE5y_bg(2,i+1)=T_cf_DICE5y_bg(2,i)+0.025*(T_cf_DICE5y_bg(1,i)-T_cf_DICE5y_bg(2,i));
end
T_cf_DICE_diff=interp1(0:5:200,T_cf_DICE5y_pulse(1,:),year) - interp1(0:5:200,T_cf_DICE5y_bg(1,:),year);

%DICE2013
T_cf_DICE2013_pulse=[[T0atm;T0ocean_DICE] zeros(2,40)];
Fpulse=FCO22x*log((M2010+Pulse)/Mequilibrium)/log(2);
for i=1:40
    T_cf_DICE2013_pulse(1,i+1)=T_cf_DICE2013_pulse(1,i)+0.098*(Fpulse-1.310*T_cf_DICE2013_pulse(1,i)-0.088*(T_cf_DICE2013_pulse(1,i)-T_cf_DICE2013_pulse(2,i)));
    T_cf_DICE2013_pulse(2,i+1)=T_cf_DICE2013_pulse(2,i)+0.025*(T_cf_DICE2013_pulse(1,i)-T_cf_DICE2013_pulse(2,i));
end
T_cf_DICE2013_pulse=transpose(interp1([0:5:200],transpose(T_cf_DICE2013_pulse),year));
%BAckground DICE2013 trajectory
T_cf_DICE2013_bg=[[T0atm;T0ocean_DICE] zeros(2,40)];
Fpulse=FCO22x*log(M2010/Mequilibrium)/log(2);
for i=1:40
    T_cf_DICE2013_bg(1,i+1)=T_cf_DICE2013_bg(1,i)+0.098*(Fpulse-1.310*T_cf_DICE2013_bg(1,i)-0.088*(T_cf_DICE2013_bg(1,i)-T_cf_DICE2013_bg(2,i)));
    T_cf_DICE2013_bg(2,i+1)=T_cf_DICE2013_bg(2,i)+0.025*(T_cf_DICE2013_bg(1,i)-T_cf_DICE2013_bg(2,i));
end
T_cf_DICE2013_bg=transpose(interp1([0:5:200],transpose(T_cf_DICE2013_bg),year));
T_cf_DICE2013_diff=T_cf_DICE2013_pulse-T_cf_DICE2013_bg;

%Identical approach for Lemoine and Rudik 
Tsteadystate=Climatesensitivity*log((M2010+Pulse)/M2010)/log(2); 
%Temperature increase=CS[log(436/289)-log(389/289)]/ln2=CS[log(436/389)]/ln2=0.51°C
%Pulse must be of the size of current excess CO2 (=M2010) in order to add 3.1°C on top of current warming (idem in Stata)
% for i=19:21
% Tdot=  @(t,T) [Tsteadystate_LR*G(i,3)/G(i,1) - G(i,3)/G(i,1)*T(1)- G(i,4)/G(i,1)*(T(1)-T(2)); ...
%     G(i,4)/G(i,2)*(T(1)-T(2))];
% twobc = @(Ta,Tb) [Ta(1) ; Ta(2)];
% solinit=bvpinit(linspace(0,600),[0.2;0.2]);
% sol = bvp4c(Tdot,twobc,solinit); %solve with ode45 would do as well since all boundary conditions are initial conditions 
% Ttemporary=deval(sol,year);
% T_cf_LR_diff(i-18,:)=Ttemporary(1,:);
% end

%----------------------------------------------------------------------------------------------------------------
%TEMPERATURE IMPACT RESPONS FUNCTION (shown in Fig 1)
T_JoosGeof_pulse=zeros(16*16,201);
for k=1:16 %loop runs over carbon decay models
  for i=1:16 %loop runs over thermal inertia models
    Forcing=Climatesensitivity*G(i,3)*log((Mequilibrium+M_JoosLRGolosGL_pulse(k,:))/Mequilibrium)/log(2);
    [t,T]=ode45(@(t,T) odefcn_inertia(t,T,Forcing,year,G(i,:)),[0 200],[T0atm;T0ocean],options);
    Ttemporary=interp1(t,T(:,1),year);
    T_JoosGeof_pulse(16*(k-1)+i,:)=transpose(Ttemporary);
  end 
end
%background temperature
T_JoosGeof_bg=zeros(16*16,201);
for k=1:16 %loop runs over carbon decay models
  for i=1:16 %loop runs over thermal inertia models
    Forcing=Climatesensitivity*G(i,3)*log((Mequilibrium+M_JoosLRGolosGL_bg(k,:))/Mequilibrium)/log(2);
    [t,T]=ode45(@(t,T) odefcn_inertia(t,T,Forcing,year,G(i,:)),[0 200],[T0atm;T0ocean],options);
    Ttemporary=interp1(t,T(:,1),year);
    T_JoosGeof_bg(16*(k-1)+i,:)=transpose(Ttemporary);
  end 
end
T_JoosGeof_diff=T_JoosGeof_pulse - T_JoosGeof_bg;
T_JoosGeof_sorted=sort(T_JoosGeof_diff,1);
T_deciles=zeros(9,201);
for i=1:9
   T_deciles(i,:)=T_JoosGeof_sorted(round(i*25.6),:);
end

%JGbestfit, LR, GHTK, GL trajectory
options= odeset('RelTol',1e-9);
for i=18:21 
Forcing=Climatesensitivity*G(i,3)*log((Mequilibrium+M_JoosLRGolosGL_pulse(i,:))/Mequilibrium)/log(2);
[t,T]=ode45(@(t,T) odefcn_inertia(t,T,Forcing,year,G(i,:)),[0 200],[T0atm;T0ocean],options);
Ttemporary=interp1(t,T(:,1),year, 'spline');
%background
Forcing=Climatesensitivity*G(i,3)*log((Mequilibrium+M_JoosLRGolosGL_bg(i,:))/Mequilibrium)/log(2);
[t,T]=ode45(@(t,T) odefcn_inertia(t,T,Forcing,year,G(i,:)),[0 200],[T0atm;T0ocean],options);
Ttemporary_bg=interp1(t,T(:,1),year, 'spline');
if i==18; T_JGbestfit_diff=Ttemporary-Ttemporary_bg;end
if i==19; T_LR_diff=Ttemporary-Ttemporary_bg;end
if i==20; T_Golosov=Ttemporary-Ttemporary_bg;
T_Golosov=interp1(3:201,T_Golosov(1,3:201), year,'spline','extrap'); % pchip replaces cubic
end
if i==21; T_GL_diff=Ttemporary-Ttemporary_bg;end
end

% Old code with the same result
% %Lemoine and Rudik 
% %pulse and dg are meaningless because they have the wrong initial conditions
% Forcing=Climatesensitivity*G(19,3)*log((M2010+M_JoosLRGolosGL_diff(19,:))/M2010)/log(2);%Pulse must be of the size of current excess CO2 (=M2015) in order to add 3.1°C
% [t,T]=ode45(@(t,T) odefcn_inertia(t,T,Forcing,year,G(19,:)),[0 200],[0;0],options);
% Ttemporary=interp1(t,T(:,1),year);
% T_LR_diff2=Ttemporary;
% sum(T_LR_diff2'-T_LR_diff')
% plot(year,T_LR_diff2,year,T_LR_diff);legend('old','new')
% %GL18
% Forcing=Climatesensitivity*G(21,3)*log((M2010+M_JoosLRGolosGL_diff(21,:))/M2010)/log(2);
% [t,T]=ode45(@(t,T) odefcn_inertia(t,T,Forcing,year,G(21,:)),[0 200],[0;0],options);
% Ttemporary=interp1(t,T(:,1),year);
% T_GL_diff2=Ttemporary;
% plot(year,T_GL_diff2,year,T_GL_diff);legend('old','new')
% %Golosov 
% T_Golosov2=Climatesensitivity*log((M2010+M_JoosLRGolosGL_diff(20,:))/M2010)/log(2);
% T_Golosov2(1,1)=T_Golosov2(1,2);  % first value is wrong because modelled with an extremely fast decay
% plot(year,T_Golosov2,year,T_Golosov);legend('old','new')

%DICE
T_DICE_pulse=[[T0atm;T0ocean_DICE] zeros(2,40)];
Fpulse=FCO22x*log(M_DICE_GtC_pulse(1,:)/588)/log(2);
for i=1:40
    T_DICE_pulse(1,i+1)=T_DICE_pulse(1,i)+0.1005*(Fpulse(1,i+1)-1.1875*T_DICE_pulse(1,i)-0.088*(T_DICE_pulse(1,i)-T_DICE_pulse(2,i)));
    T_DICE_pulse(2,i+1)=T_DICE_pulse(2,i)+0.025*(T_DICE_pulse(1,i)-T_DICE_pulse(2,i));
end
T_DICE_pulse=transpose(interp1([0:5:200],transpose(T_DICE_pulse),year));
%BAckground DICE trajectory
T_DICE_bg=[[T0atm;T0ocean_DICE] zeros(2,40)];
Fpulse=FCO22x*log(M_DICE_GtC_bg(1,:)/588)/log(2);
for i=1:40
    T_DICE_bg(1,i+1)=T_DICE_bg(1,i)+0.1005*(Fpulse(1,i+1)-1.1875*T_DICE_bg(1,i)-0.088*(T_DICE_bg(1,i)-T_DICE_bg(2,i)));
    T_DICE_bg(2,i+1)=T_DICE_bg(2,i)+0.025*(T_DICE_bg(1,i)-T_DICE_bg(2,i));
end
T_DICE_bg=transpose(interp1([0:5:200],transpose(T_DICE_bg),year));
T_DICE_diff=T_DICE_pulse-T_DICE_bg;


%DICE2100 trajectory  
scenario=3; %RCP SSP2-RCP4.5
T_DICE2100_pulse=[X2100(scenario,5:6)' zeros(2,40)];%change here number 3 for another scenario 
M2100_GtC=sum(X2100(scenario,1:4))/3.664; % above preindustrial
Fpulse=FCO22x*log((M_DICE_GtC_pulse(1,:)+M2100_GtC)/588)/log(2);
for i=1:40
    T_DICE2100_pulse(1,i+1)=T_DICE2100_pulse(1,i)+0.1005*(Fpulse(1,i+1)-1.1875*T_DICE2100_pulse(1,i)-0.088*(T_DICE2100_pulse(1,i)-T_DICE2100_pulse(2,i)));
    T_DICE2100_pulse(2,i+1)=T_DICE2100_pulse(2,i)+0.025*(T_DICE2100_pulse(1,i)-T_DICE2100_pulse(2,i));
end
T_DICE2100_pulse=transpose(interp1([0:5:200],transpose(T_DICE2100_pulse),year));
%BAckground DICE2100 trajectory
T_DICE2100_bg=[X2100(scenario,5:6)' zeros(2,40)];
Fpulse=FCO22x*log((M_DICE_GtC_bg(1,:)+M2100_GtC)/588)/log(2);
for i=1:40
    T_DICE2100_bg(1,i+1)=T_DICE2100_bg(1,i)+0.1005*(Fpulse(1,i+1)-1.1875*T_DICE2100_bg(1,i)-0.088*(T_DICE2100_bg(1,i)-T_DICE2100_bg(2,i)));
    T_DICE2100_bg(2,i+1)=T_DICE2100_bg(2,i)+0.025*(T_DICE2100_bg(1,i)-T_DICE2100_bg(2,i));
end
T_DICE2100_bg=transpose(interp1([0:5:200],transpose(T_DICE2100_bg),year));
T_DICE2100_diff=T_DICE2100_pulse-T_DICE2100_bg;

%DICE2013
T_DICE2013_pulse=[[T0atm;T0ocean_DICE] zeros(2,40)];
Fpulse=FCO22x*log(M_DICE2013_GtC_pulse(1,:)/588)/log(2); % DIce2013 had FCO22x of 3.8 but I standardise to 3.1climate sensitivity
%add a second trajectory where I add 834ppm*366GtC/47ppm in fraction above to obtain a pulse under 2.6°C  in 2100 (SSP2-RCP4.5) conditions
for i=1:40
    T_DICE2013_pulse(1,i+1)=T_DICE2013_pulse(1,i)+0.098*(Fpulse(1,i+1)-1.310*T_DICE2013_pulse(1,i)-0.088*(T_DICE2013_pulse(1,i)-T_DICE2013_pulse(2,i)));
    T_DICE2013_pulse(2,i+1)=T_DICE2013_pulse(2,i)+0.025*(T_DICE2013_pulse(1,i)-T_DICE2013_pulse(2,i));
end
T_DICE2013_pulse=transpose(interp1([0:5:200],transpose(T_DICE2013_pulse),year));
%BAckground DICE2013 trajectory
T_DICE2013_bg=[[T0atm;T0ocean_DICE] zeros(2,40)];
Fpulse=FCO22x*log(M_DICE2013_GtC_bg(1,:)/588)/log(2);
for i=1:40
    T_DICE2013_bg(1,i+1)=T_DICE2013_bg(1,i)+0.098*(Fpulse(1,i+1)-1.310*T_DICE2013_bg(1,i)-0.088*(T_DICE2013_bg(1,i)-T_DICE2013_bg(2,i)));
    T_DICE2013_bg(2,i+1)=T_DICE2013_bg(2,i)+0.025*(T_DICE2013_bg(1,i)-T_DICE2013_bg(2,i));
end
T_DICE2013_bg=transpose(interp1([0:5:200],transpose(T_DICE2013_bg),year));
T_DICE2013_diff=T_DICE2013_pulse-T_DICE2013_bg;
%Alternative calculation Gerlagh and Liski IRF 
% %identical result for decay
% a_GL=[0.163 0.184 0.449];eta_GL=[0 0.074 0.470];eps_GL=0.183;decade=year/10;
% Pi_degree100GtC=0.0156*3.664*3/0.026*0.1;% to convert %GDP damage to °C.   1.56%GDP/TtCO2 * 3.664TtCO2/TtC * 3°C/2.6%GDP * 0.1 TtC/100GtC
% M_GerlaghLiski_diff_ppm=46.9*a_GL*((1-eta_GL'*ones(1,201)).^(ones(3,1)*decade));
% %damage impact response function (bad approximation)
% T_GerlaghLiski= ones(1,3)* Pi_degree100GtC*eps_GL*...
%     [a_GL(1)*((1-eta_GL(1)).^decade-(1-eps_GL).^decade)/(eps_GL-eta_GL(1));...
%      a_GL(2)*((1-eta_GL(2)).^decade-(1-eps_GL).^decade)/(eps_GL-eta_GL(2));...
%      a_GL(3)*((1-eta_GL(3)).^decade-(1-eps_GL).^decade)/(eps_GL-eta_GL(3))];
% %Temperature IRF which gives a reasonable approximation  
%     theta_GerlaghLiski_GtC= ones(1,3)*eps_GL*... % unlike theta in their paper, pi is left out of this calculation
%     [a_GL(1)*((1-eta_GL(1)).^decade-(1-eps_GL).^decade)/(eps_GL-eta_GL(1));...
%      a_GL(2)*((1-eta_GL(2)).^decade-(1-eps_GL).^decade)/(eps_GL-eta_GL(2));...
%      a_GL(3)*((1-eta_GL(3)).^decade-(1-eps_GL).^decade)/(eps_GL-eta_GL(3))]*100; %last factor because initial function is for TC02/TtCO2injected, whereas we want GtCO2/100GtC02injected 
%     T_GerlaghLiski_bis=3.1*log(1+theta_GerlaghLiski_GtC/851)/log(2);
% clear a_GL eta_GL eps_GL decade Pi_degree100TtC
[PeakIRFtemp,I]=max([T_JGbestfit_diff ;T_FAIR_diff(1,:);T_DICE_diff(1,:);T_DICE2013_diff(1,:); T_FUND_Beta2600;T_PAGE;T_GL_diff;T_Golosov;T_LR_diff(1,:)]');
PeakIRFyear=I-1;
%PLOTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fig 1 IRF DICE PAGE FUND Lemoine&Rudik Golosov Gerlagh&Liski
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
%     p2=plot(year, T_FAIR_diff(1,:),'-b','linewidth',2);
%     %p3=plot(year,ones(1,201)*0.0017*Pulse/3.664,'--k','linewidth',1);
%     p3=plot(year, T_FAIR_diff(4,:),'-.b','linewidth',2);        
        plot(year, T_DICE_diff(1,:),'k','linewidth',2);
    p4=scatter(year(1:10:end), T_DICE_diff(1,1:10:end),'s','MarkerEdgeColor',[0.25 0.6 1],'linewidth',1.5);
        plot(year, T_DICE2013_diff(1,:),'k','linewidth',1.5);
    p5=scatter(year(1,1:10:end), T_DICE2013_diff(1,1:10:end),'p','MarkerEdgeColor','c','linewidth',1.5);
        plot(year, T_FUND_Beta2600,'k','linewidth',1.5);
    p6=scatter(year(1:10:end), T_FUND_Beta2600(1,1:10:end),'v','MarkerEdgeColor','m','linewidth',1.5); 
        plot(year, T_PAGE,'k','linewidth',1.5); 
    p7=scatter(year(1:10:end), T_PAGE(1:10:end),'d','MarkerEdgeColor','r','linewidth',1.5);  
         plot(year, T_Golosov, 'k','linewidth',1.5);
    p8=scatter(year(1:10:end), T_Golosov(1:10:end), '+','MarkerEdgeColor',[0.9 0.55 0],'linewidth',1.5);
        plot(year, T_GL_diff, 'k', 'linewidth', 1.5);        
    p9=scatter(year(1:10:end), T_GL_diff(1:10:end), 'x','MarkerEdgeColor',[1 0.9 0.1],'linewidth', 1.5);
        plot(year, T_LR_diff(1,:),'k','linewidth',1.5);
    p10=scatter(year(1:10:end), T_LR_diff(1,1:10:end),'^','MarkerEdgeColor','b','linewidth',1.5);
    hold off
    legend([ar p1 p4 p5 p6 p7 p8 p9 p10],'Deciles CMIP5 combinations','Best fit CMIP5 ensemble',...
        'DICE16','DICE13','FUND','PAGE','GHKT14','GL18','LR17',...
    'Location','SouthEast','FontSize',11,'NumColumns',2) %,'Orientation','horizontal' or 'Southoutside'
    axis([0 200 0 0.0032*Pulse/3.664])  %0.32 for 100GtC pulse
    xlabel('years')
    ylabel('Temperature increase (°C)')
    pbaspect([1.5 1 1])
    set(gca,'fontsize',11)
    set(gcf,'PaperOrientation','landscape');
    print('-fillpage','Fig1', '-dpdf'); % '-fillpage' gives a too high figure
%----------------------------------------------------------------------------------------------------------------------------------
%Fig 4 Thermal inertia
    subplot(1,1,1)
%     plot(year,T_cf_JoosGeof_LR_GL_diff(1:16,:),'color',[0 0.5 0.25],'linewidth',0.1);
    p2=plot(year,T_cf_JoosGeof_LR_GL_diff(18,:),'r','linewidth',3);
    hold on
    plot(year,T_cf_DICE_diff,'k','linewidth',1);
    p3=scatter(year(1,1:10:end),T_cf_DICE_diff(1,1:10:end),'s','MarkerEdgeColor',[0.25 0.6 1],'linewidth',1.5);
    plot(year,T_cf_DICE2013_diff(1,:),'k','linewidth',1);
    p4=scatter(year(1,1:10:end),T_cf_DICE2013_diff(1,1:10:end),'p','MarkerEdgeColor','c','linewidth',1.5);
    plot(year,T_cf_FUND_diff,'k','linewidth',1);
    p5=scatter(year(1,1:10:end),T_cf_FUND_diff(1,1:10:end),'v','MarkerEdgeColor','m','linewidth',1.5);
    plot(year,T_cf_PAGE_diff,'k','linewidth',1);
    p6=scatter(year(1,1:10:end),T_cf_PAGE_diff(1,1:10:end),'d','MarkerEdgeColor','r','linewidth',1.5); 
    plot(year(1:10:end),ones(1,21)*Tsteadystate,'k','linewidth',1); %Golosov et al.
    p7=scatter(year(1:10:end),ones(1,21)*Tsteadystate,'+','MarkerEdgeColor',[0.9 0.55 0],'linewidth',1.5);  %Golosov
    plot(year,T_cf_JoosGeof_LR_GL_diff(21,:),'k', 'linewidth', 1);
    p8=scatter(year(1,1:10:end),T_cf_JoosGeof_LR_GL_diff(21,1:10:end),'x','MarkerEdgeColor',[0.25 0.6 1], 'linewidth', 1.5);
    plot(year,T_cf_JoosGeof_LR_GL_diff(19,:),'k','linewidth',1); %LR
    p9=scatter(year(1,1:10:end),T_cf_JoosGeof_LR_GL_diff(19,1:10:end),'^','MarkerEdgeColor','b','linewidth',1.5); %LR
    hold off
    %title('Increase in temperature for a constant increase in CO2 concentration by 47 ppm')
    axis([0 200 0 0.55])
    legend([ p2; p3; p4; p5; p6; p7; p8;p9],'Best fit CMIP5 ensemble','DICE16','DICE13','FUND','PAGE','GHKT14','GL18',...
        'LR17','Location','SouthEast', 'FontSize',11,'NumColumns',2);
    xlabel('years')
    ylabel('Temperature increase (°C)')
    set(gca,'fontsize',11)
    set(gcf,'PaperOrientation','landscape');
    print('-fillpage','Fig4', '-dpdf');
%Fig 3 Decay All models
    subplot(1,1,1);
    %plot(year,M_JoosLRGolosGL_diff(1:16,:)*46.9/366.3,'color',[0 0.5 0.25],'linewidth',0.1);
    p2=plot(year,M_JoosLRGolosGL_diff(18,:)*46.9/366.3,'r','linewidth',3);
    hold on
    p3=plot(year,M_FAIR_diff(1,:)*46.9/366.3,'--r','linewidth',1);
        %p3=scatter(year(1,1:10:end),M_FAIR_diff(1,1:10:end)*46.9/366.3,'o','linewidth',1);
    p4=plot(year,M_FAIR_diff(4,:)*46.9/366.3,':r','linewidth',1);
        %p4=scatter(year(1,1:10:end),M_FAIR_diff(4,1:10:end)*46.9/366.3,'o','linewidth',1);
    plot(year,M_DICE_diff*46.9/366.3,'k','linewidth',1);
        p5=scatter(year(1,1:10:end),M_DICE_diff(1,1:10:end)*46.9/366.3,'s','MarkerEdgeColor',[0.25 0.6 1],'linewidth',1.5);
    plot(year(1,1:5:end),M_DICE2013_diff(1,1:5:end)*46.9/366.3,'k','linewidth',1);
        p6=scatter(year(1,1:10:end),M_DICE2013_diff(1,1:10:end)*46.9/366.3,'p','MarkerEdgeColor','c','linewidth',1.5);
    plot(year,M_FUND_diff_ppm,'k','linewidth',1);
        p7=scatter(year(1,1:10:end),M_FUND_diff_ppm(1,1:10:end),'v','MarkerEdgeColor','m','linewidth',1.5);
    plot(year,M_PAGE_diff_ppm,'k','linewidth',1);  
        p8=scatter(year(1,1:10:end),M_PAGE_diff_ppm(1,1:10:end),'d','MarkerEdgeColor','r','linewidth',1.5);
    plot(year(1:10:end),M_JoosLRGolosGL_diff(20,1:10:end)*46.9/366.3,'k','linewidth',1);
        p9=scatter(year(1:10:end),M_JoosLRGolosGL_diff(20,1:10:end)*46.9/366.3,'+','MarkerEdgeColor',[0.9 0.55 0],'linewidth',1.5);
     plot(year,M_JoosLRGolosGL_diff(21,:)*46.9/366.6,'k', 'linewidth', 1);
        p10=scatter(year(1,1:10:end),M_JoosLRGolosGL_diff(21,1:10:end)*46.9/366.6,'x','MarkerEdgeColor',[01 0.9 0.1], 'linewidth', 1.5);
    plot(year,M_JoosLRGolosGL_diff(19,:)*46.9/366.3,'k','linewidth',1);
        p11=scatter(year(1,1:10:end),M_JoosLRGolosGL_diff(19,1:10:end)*46.9/366.3,'^','MarkerEdgeColor','b','linewidth',1.5);
    hold off
    axis([0 200 0 50])
    legend([ p2 p3 p4 p5 p6 p7 p8 p9 p10 p11],'Best fit CMIP5 ensemble','FAIR','FAIR in 2100 (SSP2-RCP 4.5)',...
        'DICE16','DICE13','FUND','PAGE','GHKT14','GL18','LR17', 'FontSize',11,'NumColumns',2)
    xlabel('years')
    ylabel('CO2 concentration increase (ppm)')
    pbaspect([1.5 1 1])
    set(gca,'fontsize',11)
    set(gcf,'PaperOrientation','landscape');
    print('-fillpage','Fig3', '-dpdf');
%IRF 2010 vs 2100. Not reported in the paper.
    subplot(1,1,1);
    ar=area(transpose(year),transpose(T_deciles(9,:)), 'FaceColor' , [0 0.90 0.45]);
    hold on
    area(transpose(year),transpose(T_deciles(8,:)), 'FaceColor' , [0 0.7 0.35])
    area(transpose(year),transpose(T_deciles(7,:)), 'FaceColor' , [0 0.5 0.25])
    area(transpose(year),transpose(T_deciles(6,:)), 'FaceColor' , [0 0.4 0.2])
    area(transpose(year),transpose(T_deciles(4,:)), 'FaceColor' , [0 0.5 0.25])
    area(transpose(year),transpose(T_deciles(3,:)), 'FaceColor' , [0 0.7 0.35])
    area(transpose(year),transpose(T_deciles(2,:)), 'FaceColor' , [0 0.90 0.35])
    area(transpose(year),transpose(T_deciles(1,:)), 'FaceColor' , [1 1 1])
    p1=plot(year,T_JGbestfit_diff,'color',[0 0.2 0],'linewidth',1);
    p2=plot(year, T_FAIR_diff(1,:),'b','linewidth',1);
    p3=plot(year, T_DICE_diff(1,:),'r','linewidth',1);
    p4=plot(year, T_DICE2100_diff(1,:),'-.r','linewidth',1);
    %p4=plot(year, T_FAIR_diff(2,:),'--r','linewidth',1);
    %p5=plot(year, T_FAIR_diff(3,:),'--r','linewidth',1);
    p6=plot(year, T_FAIR_diff(4,:),'-.b','linewidth',1);
    %p7=plot(year, T_FAIR_diff(5,:),':r','linewidth',1);
    hold off
    title('Increase in temperature after a 100GtC pulse of emissions')
    legend([ar p1 p2 p6 p3 p4],'Deciles CMIP5 combinations','Best fit CMIP5 ensemble','FAIR in 2010','FAIR in 2100 (SSP2-RCP 4.5)' ,'DICE in 2010',...
        'DICE in 2100 (SSP2-RCP 4.5)','Location','SouthEast')
    axis([0 200 0 0.32])
    xlabel('years')
    ylabel('Temperature increase (°C)')
    set(gcf,'PaperOrientation','landscape');
    print('-fillpage','IRF100GtC_2010_2100', '-dpdf');
%------------------------------------------------------------------------------------------
%Graphs in the AER comment (2020) to Lemoine and Rudik AER (2017)
subplot(3,1,1);
    plot(year,transpose(M_JoosLRGolosGL_diff(1:16,:)*46.9/366.3),'color',[0 0.5 0.25]);
    hold on
    p1=plot(year,transpose(M_JoosLRGolosGL_diff(1,:)*46.9/366.3),'color',[0 0.5 0.25]);
    p2=plot(year,transpose(M_JoosLRGolosGL_diff(18,:)*46.9/366.3),'r','linewidth',1);
    p3=plot(year,transpose(M_JoosLRGolosGL_diff(19,:)*46.9/366.3),'b','linewidth',1);
    hold off
    title('Increase in atmospheric CO2 concentrations after a 100GtC pulse of emissions (initial increase of 47 ppm)')
    axis([0 200 0 inf])
    legend([p1 p2 p3],'CMIP5 ensemble (Joos et al.)','Best fit CMIP5 ensemble','Lemoine & Rudik')
    xlabel('years')
    ylabel('CO2 concentration increase (ppm)')
subplot(3,1,2);
    plot(year,transpose(T_cf_JoosGeof_LR_GL_diff(1:16,:)),'color',[0 0.5 0.25]);
    hold on
    p1=plot(year,transpose(T_cf_JoosGeof_LR_GL_diff(1,:)),'color',[0 0.5 0.25]);
    p2=plot(year,transpose(T_cf_JoosGeof_LR_GL_diff(17,:)),'r','linewidth',1); 
    p3=plot(year,T_cf_JoosGeof_LR_GL_diff(19,:),'b','linewidth',1);
    p4=plot(year,T_cf_JoosGeof_LR_GL_diff(20,:),year,T_cf_JoosGeof_LR_GL_diff(21,:),'Color','b','linewidth',0.5);
    hold off
    title('Increase in temperature for a constant increase in CO2 concentration by 47 ppm (for a common 3 climate sensitivity)')
    axis([0 200 0 0.5])
    legend([p1 p2 p3],'CMIP5 ensemble (Geoffroy et al.)','Best fit CMIP5 ensemble','Lemoine & Rudik', 'Location','SouthEast')
    xlabel('years')
    ylabel('Temperature increase (°C)')
subplot(3,1,3);
    ar=area(transpose(year),transpose(T_deciles(9,:)), 'FaceColor' , [0 0.90 0.45]);
    hold on
    area(transpose(year),transpose(T_deciles(8,:)), 'FaceColor' , [0 0.70 0.35])
    area(transpose(year),transpose(T_deciles(7,:)), 'FaceColor' , [0 0.5 0.25])
    area(transpose(year),transpose(T_deciles(6,:)), 'FaceColor' , [0 0.2 0.2])
    area(transpose(year),transpose(T_deciles(4,:)), 'FaceColor' , [0 0.5 0.25])
    area(transpose(year),transpose(T_deciles(3,:)), 'FaceColor' , [0 0.5 0.25])
    area(transpose(year),transpose(T_deciles(2,:)), 'FaceColor' , [0 0.70 0.35])
    area(transpose(year),transpose(T_deciles(1,:)), 'FaceColor' , [1 1 1])
    p1=plot(year,T_JGbestfit_diff,'r','linewidth',1);
    p2=plot(year, T_LR_diff(1,:),'b','linewidth',1);
    plot(year, T_LR_diff(2,:),year, T_LR_diff(3,:),'Color','b','LineWidth',0.5);
    hold off
    title('Increase in temperature after a 100GtC pulse of emissions')
    legend([ar p1 p2],'Deciles CMIP5 combinations','Best fit CMIP5 ensemble','Lemoine & Rudik','Location','SouthEast')
    xlabel('years')
    ylabel('Temperature increase (°C)')
    print('-fillpage','IRF100GtC_2010_LR', '-dpdf'); % study works wit climate sensitivity of 3°C, change above
