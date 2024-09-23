clc;
clear;
close all;
%% define paths
SetC.interp=1; %interpolation needed
%scebario 1=2020, 2=2025, 3 =2030
Scenario=1;
EV=1;
Public_flag=1;% Public charging station flag
if Scenario==1
    SetC.netcase= 'MG';
end
path=Allp(SetC);
Scenario=1;
addpath(genpath(path.dat1));%Path for the files and folder
%% Load Data
run LoadDat;
%% PARAMET
yalmip('clear');
opt.Horizon   =8760;
opt.Period    = 10; %number of years
opt.scenarios = 1; %number of years
%% Interpolate Profile Data AND LOAD PROFILE ALLOCATION
PROF=InterpDat(PROF,ANZ,SetC);                                              %interpolate Profile Data
[NOD,PROF]=AllocProfs(GEN,LOAD,PROF,ANZ,NOD,EV);                               %allocate profile data to nodes
BRANCH=EditBranch(ANZ,NOD,LIN,TRA,SWT);
[NET,indK]=Admitt(BRANCH,ANZ,NOD,LIN,TRA,SWT);
mpc.Vbase=0.4*1e3;                                                          %V
mpc.baseMVA=0.630*1e6;                                                       %KVA to VA
mpc.Zbase=mpc.Vbase.^2 /mpc.baseMVA;
mpc.Ibase=mpc.Vbase /mpc.Zbase;
%% form index lists for slack, PV, and PQ buses
%% tariff
cost.penalty1 =0.11;                                                        %curtailed penalty
cost.choose_tarrif=2;
load('MG.mat')
if cost.choose_tarrif==2
    cost.WT= wtarrifs(NOD.PLoadProf,opt.Horizon);                                      %make it Euro/kwh
     cost.price=cost.WT(1,1:opt.Horizon)./100;  
else
    cost.fixed=0.30;                                                        %euro/kwhr
    cost.price=cost.fixed*ones(1,opt.Horizon);
end
%% wind system
%% Wind
% Specify filenames
%% Wind
load('wind.mat', 'wind');
Grid.wind=wind';
Grid.wind_capacity=50;      % capacity limit kw
Grid.wind_eff=0.18;          % pv efficiency
Grid.N_wind=1;
Grid.windIncMatrix= zeros(ANZ.K,Grid.N_wind);
wind_C = [1];
for i = 1: Grid.N_wind
Grid.windIncMatrix(wind_C (i,1),i) = 1;
end
%Ppv(:,1:opt.Horizon).*1e3/mpc.baseMVA;
% Capacity variables
Grid.wind_p=sdpvar(1,opt.Horizon);
Grid.wind_cap=sdpvar(1,opt.Period);
%Grid.PV_irr=Grid.pv_eff*Grid.pv_irradiation(1,97:120);
%Grid.PV_irr=Grid.pv_irradiation(1,1:opt.Horizon);
Grid.wind=Grid.wind(1,1:opt.Horizon).*0.5;
%% PV syastem
%% PV profile
%Ppv=sum(NOD.PGenProf);
Grid.pv_capacity=500;      % capacity limit kw
Grid.pv_eff=0.18;          % pv efficiency
% Solar irradiation in kW/m^2
Grid.pv_irradiation=PROF.PGenProf1'; % in kW/m^2
Grid.N_PV=1;
Grid.PV_PIncMatrix= zeros(ANZ.K,Grid.N_PV);
Grid.PV_P = [1];
for i = 1: Grid.N_PV
Grid.PV_PIncMatrix(Grid.PV_P(i,1),i) = 1;
end
Grid.PV=sdpvar(1,opt.Horizon);
%Ppv(:,1:opt.Horizon).*1e3/mpc.baseMVA;
% Capacity variables
Grid.pv_cap=sdpvar(1,opt.Period);
%Grid.PV_irr=Grid.pv_eff*Grid.pv_irradiation(1,97:120);
Grid.PV_irr=Grid.pv_irradiation(1,1:opt.Horizon);
%% grid model
Grid.Pload  = ((NOD.PLoadProf(:,1:opt.Horizon)));%.*1e3./mpc.baseMVA;         %active load in per unit in watts
Grid.Qload  = (NOD.QLoadProf(:,1:opt.Horizon).*1e3)./mpc.baseMVA;           %reactive load  in per unit in watts


Grid.V      = sdpvar(ANZ.K,opt.Horizon);                            %node voltages^2 
Grid.S      = sdpvar(ANZ.L,opt.Horizon);
Grid.R_Line      = (LIN.l.*LIN.r)./mpc.Zbase;
Grid.X_Line      = (LIN.l.*LIN.x)./mpc.Zbase;
Grid.L      = sdpvar(ANZ.L,opt.Horizon);
Grid.P_in=sdpvar(ANZ.K,opt.Horizon);
Grid.Pshed= sdpvar(1,opt.Horizon);% grid supplied P
Grid.Pshed_cap= sdpvar(1,opt.Period);% grid supplied P
Grid.Pgmax=[ones(1,opt.Horizon);zeros(ANZ.K-1,opt.Horizon)];
Grid.Qgmax=[ones(1,opt.Horizon);zeros(ANZ.K-1,opt.Horizon)];
%initial time variables for 24hrs
dt=1; % Time step
Grid.inc=full(sparse(NET.vK',1:size(NET.vK,1),1,ANZ.K,size(NET.vK,1)));
Grid.finc=full(sparse(NET.zK',1:size(NET.zK,1),1,ANZ.K,size(NET.zK,1)));
Grid.Pshed_cap=sdpvar(1,opt.Period);
%% making ambiguity set for EV based load

%%  Battery 
N_bat=1;
battSocMin = 0;
battSocMax = 1;
%Bat.battPmax = sdpvar(N_bat,opt.Period,'full');
Bat.bat_decay=0.001;                                                            % decay rate (storage losses) %/h/100
Bat.bat_loss=0.1;                                                               % charging/discharging loss rate %/100
Bat.battCapa =sdpvar(N_bat,opt.Period);                                           % [Wh] - from sizing
%assign(Bat.battCapa,50);
battEff = 0.9;
Bat.battChargeOrNot = binvar(2,opt.Horizon);                                 % binary variable for charge or discharge
Bat.Pbatt_ch = sdpvar(N_bat,opt.Horizon);  
Bat.Pbatt_dis = sdpvar(N_bat,opt.Horizon);  
%Bat.Pbatt_pch = sdpvar(N_bat,opt.Period,opt.Horizon,'full'); 
Bat.bat_store = sdpvar(N_bat,opt.Horizon);
%soc = sdpvar(1,Horizon,'full'); 
Bat.bat_store(N_bat,1)=0;
%Bat.bat_capital1=100;
dt = 1;
Bat.batIncMatrix = zeros(ANZ.K,N_bat);
bat_gen = [1];
%index PowersystemBus
for i = 1: N_bat
    Bat.batIncMatrix(bat_gen(i,1),i) = 1;
end

%% heating netwrok modelling
%% heating load 
load(path.space_heating);
load(path.water_heating);
heat.heat_load=((space_heat_minute).*4);%+water_heat_minute); % kW
PROF.Hload=InterpDat_heat(heat.heat_load,ANZ,SetC);% interpolate Profile Data
heat.Hload=PROF.Hload';
heat.Hload=heat.Hload(:,1:opt.Horizon);
%%  thermal 
N_th=1;
heat.th_decay=0.01;                 % decay rate (storage losses) %/h/100
heat.th_loss=0.1;                   % charging/discharging loss rate %/100
heat.th_Capa =sdpvar(N_th,opt.Period);   %224; % [Wh] - from sizing
%assign(heat.th_Capa,1000);
heat.th_in=sdpvar(N_th,opt.Horizon);
heat.th_out=sdpvar(N_th,opt.Horizon);
heat.th_store=sdpvar(N_th,opt.Horizon);
heat.thChargeOrNot=binvar(2,opt.Horizon); % binary variable for charge or discharge
heat.th_store(N_th,1)=0;
heat.dt = 1;

%% heat network
n_T=opt.Horizon;
heat.HeatBranch=xlsread(path.heat_xlsx,'Branch');
heat.HeatBus=xlsread(path.heat_xlsx,'Bus');
heat.HeatBranch(:,4) =heat.HeatBranch(:,4);%(kg/h)
heat.n_HeatBranch = size(heat.HeatBranch,1);
heat.n_HeatBus = size(heat.HeatBus,1);
%heat.HeatFlowInMatrix = zeros(n_HeatBus,n_HeatBranch);  % With flow
%heat.HeatFlowInIncMatrix = zeros(n_HeatBranch,n_HeatBus);% Incidence matrix
heat.TmprtrFromDir = sdpvar(heat.n_HeatBranch, n_T);  
heat.TmprtrToDir = sdpvar(heat.n_HeatBranch, n_T);   
heat.TmprtrFromRev = sdpvar(heat.n_HeatBranch, n_T);  
heat.TmprtrToRev = sdpvar(heat.n_HeatBranch, n_T);    
heat.TmprtrBusDir = sdpvar(heat.n_HeatBus,n_T);      
heat.TmprtrBusRev = sdpvar(heat.n_HeatBus,n_T);      
%heat.FromDir = sdpvar(n_HeatBranch, Horizon);  %Temperature of branch head junction in positive direction
%heat.ToDir = sdpvar(n_HeatBranch, Horizon);    %End junction temperature of branch in positive direction
%heat.FromRev = sdpvar(n_HeatBranch, Horizon);  %Reverse branch junction temperature
%heat.ToRev = sdpvar(n_HeatBranch, Horizon);    %Temperature at the end of the reverse branch
%heat.BusDir = sdpvar(n_HeatBus,Horizon);       %The temperature of hot water at each node of the system in the positive direction
%heat.BusRev = sdpvar(n_HeatBus,Horizon);       %The temperature of hot water at each node of the reverse system
%heat.HeatCHP = sdpvar(n_CHPgen,Horizon);             %chp unit heat outpu
%heat.HeatEBoiler = sdpvar(n_EBoiler,Horizon);        %electric boiler heat output
%heat.PowerEBoiler = sdpvar(n_EBoiler,Horizon);       %Electric boiler power consumption
%heat.HeatD= sdpvar(n_HeatBus,Horizon);
heat.Cp = 4200;% joule/kgC
heat.Cp=heat.Cp/3600000;% change into Kwhr/kgC   %%%%%%%%%/3600000; % j into wh/kg*C.....4200;%/3600000/1000;  %J/(kg*C¡ã)->MW*h/(kg*C¡ã)
heat.SituationTempreture = [
    -10 -10 -8.84   -9.42   -9.42   -9.42   -8.84   -8.26   -7.10   -6.52   -5.94   -5.35   -4.77   -4.77   -4.77   -5.35   -5.94   -6.52   -6.52   -6.52   -7.10   -7.68   -8.26   -8.26
];
%SituationTempreture=repelem(SituationTempreture,1,1460);
heat.HeatD=sdpvar(heat.n_HeatBus,opt.Horizon);
%% heat pump
heat.N_HP=1;
heat.HeatSource = sdpvar(heat.n_HeatBus,opt.Horizon);        %Heat source heating
heat.hp_in=sdpvar(heat.N_HP,opt.Horizon);
heat.hp_out=sdpvar(heat.N_HP,opt.Horizon);
heat.hp_eff=3.2;                          % power to heat COP
heat.HP_life_time=20;
heat.hp_capacity=10;                    % capacity limit kW

heat.hp_cap=sdpvar(heat.N_HP,opt.Period);
%assign(heat.hp_cap,500);

%% HeatBus PowersystemBus HP0
%bus(Heat)  bus(power system)
heat.powerHeatpumpIncMatrix= zeros(ANZ.K,heat.N_HP);
heat.SourceHeatpumpIncMatrix=zeros(heat.n_HeatBus,heat.N_HP);
heat.E_HP = [1 3];
for i = 1: heat.N_HP
    heat.SourceHeatpumpIncMatrix(heat.E_HP(i,1),i) = 1;
    heat.powerHeatpumpIncMatrix(heat.E_HP(i,2),i)=1;
end
heat.thIncMatrix = zeros(heat.n_HeatBus,N_th);
heat.thoutMatrix = zeros(heat.n_HeatBus,N_th);
heat.th_gen = [3 1];
%index heatystemBus
for i = 1: N_th
    heat.thIncMatrix(heat.th_gen(i,1),i) = 1;
    heat.thoutMatrix(heat.th_gen(i,2),i)=1;
end
heat.Heat_load=zeros(heat.n_HeatBus,opt.Horizon);
heat.Heat_load((heat.HeatBus(:,2)==2),:)=heat.Hload;
heat.HeatFlowInMatrix = zeros(heat.n_HeatBus,heat.n_HeatBranch);       
heat.HeatFlowInIncMatrix = zeros(heat.n_HeatBranch,heat.n_HeatBus);    
for i=1:heat.n_HeatBranch
    %Tobus
    heat.HeatFlowInMatrix(heat.HeatBranch(i,3),i) = 1*heat.HeatBranch(i,4);
    heat.HeatFlowInIncMatrix(i,heat.HeatBranch(i,3)) = 1;
end
heat.HeatFlowInBus = heat.HeatFlowInMatrix*ones(heat.n_HeatBranch,1);   
heat.HeatFlowOutMatrix = zeros(heat.n_HeatBus,heat.n_HeatBranch);       
heat.HeatFlowOutIncMatrix = zeros(heat.n_HeatBranch,heat.n_HeatBus);    
for i=1:heat.n_HeatBranch
    %Frombus
    heat.HeatFlowOutMatrix(heat.HeatBranch(i,2),i) = 1*heat.HeatBranch(i,4);
    heat.HeatFlowOutIncMatrix(i,heat.HeatBranch(i,2)) = 1;
end
%% mas flow is in kg per hour
heat.HeatFlowOutBus = heat.HeatFlowOutMatrix*ones(heat.n_HeatBranch,1);
%% heat line losses
    heat.coefficient = zeros(heat.n_HeatBranch,1);
    for i = 1: heat.n_HeatBranch
        %coefficient(i) = exp(-HeatBranch(i,8)*HeatBranch(i,5)/4200/HeatBranch(i,4)*3600);
       heat.coefficient(i) = exp(-heat.HeatBranch(i,8)*heat.HeatBranch(i,5)/heat.Cp/heat.HeatBranch(i,4)/100000);
    end
%==========economy parameters
a=0.06;  %interest rate, residual value and proposed lifecycle
v=1;
Ta=50;
Tloss=0.1; grideff=0.80;
AF1=a*((1+a)^Ta); af1=((1+a)^Ta)-1;
AF = AF1/af1;  %annuity factor
    %% Electrolyzer  and fuel cell
% ==============parameters variables for electrolyser============
Grid.N_elez=1;
Grid.elez=10; %specific electrical energy demand of electrolyser
elc_fc_lif_time=60000;% hours
Grid.eleyr=elc_fc_lif_time/8760; %life span of electrolyser
Grid.elzcost=300; %unit cost of elz investment
Grid.elzmc=0.08; %proportion of unit cost for maintenance
Grid.elzrep=(Ta/Grid.eleyr)-1;
Grid.Pelez=sdpvar(Grid.N_elez,opt.Horizon);
Grid.ELEZ_cap=sdpvar(1,opt.Period);

Grid.hLHV=33.3; %low heating value of hydrogen gas in kw/m3
Grid.Helz=sdpvar(Grid.N_elez,opt.Horizon);
Grid.Pemeff= 0.35; Grid.sofeff=0.54; Grid.Mofceff= 0.47; %electrical efficiency of FC
Grid.phloss=0.14; Grid.sofhloss=0.26; Grid.mofhloss=0.18; %heat loss of considered FC
Grid.pemyr=elc_fc_lif_time/8760; Grid.sofyr=20; Grid.mofyr=25; %life span of each FC
Grid.pemcost=650; Grid.sofcost=300; Grid.mofcost=350; %unit cost of investment
Grid.pemmc=0.08; Grid.sofmc=0.08; Grid.mofmc=0.08; % proportion of unit cost for maintenance
Grid.pemrep=(Ta/Grid.pemyr)-1; Grid.sofrep=(Ta/Grid.sofyr)-1; Grid.mofrep=(Ta/Grid.mofyr)-1;

Grid.h1=sdpvar(Grid.N_elez,opt.Horizon); 
Grid.Pfc=sdpvar(Grid.N_elez,opt.Horizon);
Grid.fc_cap=sdpvar(1,opt.Period);
Grid.Qfc=sdpvar(Grid.N_elez,opt.Horizon);

%% HeatBus PowersystemBus HP
%bus(Heat)  bus(power system)
Grid.powerelezIncMatrix= zeros(ANZ.K,Grid.N_elez);
Grid.SourceHeatfcIncMatrix= zeros(heat.n_HeatBus,heat.N_HP);
Grid.powerfcIncMatrix= zeros(ANZ.K,Grid.N_elez);
E__elez = [4 1 1];
for i = 1: Grid.N_elez
    Grid.powerelezIncMatrix(E__elez(i,1),i) = 1;
    Grid.powerfcIncMatrix(E__elez(i,2),i) = 1;
    Grid.SourceHeatfcIncMatrix(E__elez(i,3),i) = 1;
end
%% hydrogen storage
%============== hygrogen storage variables==========
Grid.echh2=0.98; Grid.edchh2=0.98;
Grid.h2chmax=700; % Battery max allowed charging power in kW
Grid.h2dchmax=800; % Battery max allowed discharging power in kW
Grid.selfdchh=0.01; %self-discharging parameters
cost.H2_capital=35; cost.HSmc=0.03; cost.HSyr=10; %TES unit cost, maintenance and life span
Grid.HSdgcost=0.1; %TES charging and discharging degradation cost
Grid.HSrep=(Ta/cost.HSyr)-1;
Grid.H2ChargeOrNot = binvar(2,opt.Horizon);        
cost.H2_op_cost=cost.H2_capital*0.01;
Grid.H2_store=sdpvar(N_th,opt.Horizon);     
Grid.h2ch=sdpvar(N_th,opt.Horizon);   
Grid.h2dch=sdpvar(N_th,opt.Horizon);
Grid.H2_Capa =sdpvar(N_th,opt.Period);   %224; % [Wh] - from sizing
Grid.H2_store(1,1)=0;
%% CHP
CHP.n_CHP=1;
CHP.P_CHP= sdpvar(1,opt.Horizon);
CHP.Heat_CHP= sdpvar(CHP.n_CHP,opt.Horizon);
CHP.SourceCHPheatIncMatrix = zeros(heat.n_HeatBus,CHP.n_CHP);
CHP.CHPgen = [
%HeatBus PowersystemBus
    1       4];
for i = 1: CHP.n_CHP
    CHP.SourceCHPheatIncMatrix(CHP.CHPgen(i,1),i) = 1;
end

CHP.P_CHP_cap=sdpvar(1,opt.Period);
CHP.life_time=30;

CHP.CHPIncMatrix = zeros(ANZ.K,CHP.n_CHP);
chp_gen = [1];
%index PowersystemBus
for i = 1: N_bat
    CHP.CHPIncMatrix(chp_gen(i,1),i) = 1;
end
%% start of optimization
tic;
%% monte carlo simulation for EV travelling and charging and travelling behavior
%   if EV==1
%       max_charging_energy=22*24; % maximium charging energy per day
%       EVCS.EVCS_capital=100; % 100 euros per EVCS per year cost
%       EVCS.EVCS_decision=intvar(1,opt.Period); % EVCS decision varaable
%       EVCS.PCS=22;%charging power  kw
%       EVCS.chargetime=2;%charging time
%       EVCS.nEV=100; % number of EV
%       EVCS.EVCS_max=5;
%       EV_behaviour=EV_sim(EVCS.PCS,EVCS.nEV,EVCS.chargetime);% monte carlo simulation
%       temp1=double(double(EV_behaviour.ev_charged));
%       t1=1:1:525600;
%       t2=1:60:525600;
%      for p=1:EVCS.nEV
%         EV_behaviour.EV_LP(p,:)=interp1(t1,temp1(p,:),t2).*EVCS.PCS;
%      end 
%   end
%load('EV2.mat')
% Grid.alpha=sdpvar(1,opt.Period);
% Grid.nEVCS_Uncer = sdpvar(1,opt.Period);
% Grid.EVCS_decision = sdpvar(1,opt.Period);
Grid.EVCS_max=4;
Grid.nEV=5;
Grid.EVCS_inc = [3];
Grid.delta = 0.5;
%index PowersystemBus
Grid.N_EVCS=1;
Grid.B=binvar(1,Grid.EVCS_max);
Grid.EVCSIncMatrix= zeros(ANZ.K,Grid.N_EVCS);
for i = 1:Grid.N_EVCS
    Grid.EVCSIncMatrix(Grid.EVCS_inc(i,1),i) = 1;
end
%Grid.nEVCS=[2,2,3,3,4];
%Grid.EVCS_load=sdpvar(ANZ.K,opt.Horizon,opt.Period,'full');
%% Defining variables necessary for binary linearization
%Grid.bin_max = round(log2(Grid.EVCS_max)+1);
%Grid.v = sdpvar(1,Grid.bin_max);
%Grid.z= binvar(1,Grid.bin_max);
%% Critical cost of the deterministic planning model
%Grid.co=[10,11,12,13,14];
%% EVCS 
Grid.EVCS_cap=sdpvar(1,opt.Period);
%Grid.EVCS_dec=sdpvar(1,opt.Period,'full'); 
%EVCS_load=sdpvar(1,opt.Horizon,'full');
%Grid.EVCS_dec=((sum(sum(EV_behaviour.v_driving_km_pp(1:7,:))./60))*.179).*0.15;
%Grid.EVCS_load=sum(EV_behaviour.EV_LP(1:7,:));
Private_per=0.5;
%algo=2
%% cost model
%% cost model
cost.C_gas=0.054;% euros per kwhr
cost.wind_capital=73.9;     % capital costs euros/kw
cost.wind_op_cost=cost.wind_capital*0.01;
cost.pv_capital=40;     % capital costs euros/kw
cost.pv_op_cost=cost.pv_capital*0.01;
cost.bat_capital=66;
cost.inv_cost=175/8;
cost.inv_op=cost.inv_cost*.01;
cost.bat_op_cost=cost.bat_capital*0.01;
%Wind power plans73.95 euro/kw
cost.th_capital=21.1;                 % capital costs euro/kWh
cost.th_op_cost=cost.th_capital*0.01;
cost.hp_capital=40;%4000/heat.HP_life_time;                       % capital costs Euro/kW
cost.hp_op_cost=cost.hp_capital*0.01;
cost.ELEZ_capital=Grid.elzcost.*2/Grid.eleyr;
cost.ELEZ_ope=cost.ELEZ_capital*0.01;
cost.CHP_inves_cost=1000/CHP.life_time;% euro /per kW per annual
cost.CHP_op_cost=cost.CHP_inves_cost*0.01;
cost.EVCS_capital=10000;
cost.EVCS_op_cost=cost.EVCS_capital*0.01;
cost.fc_capital=Grid.pemcost./Grid.pemyr;
cost.fc_op_cost=cost.fc_capital*0.01;
cost.waiting_cost = 20;
%% CO2 emission metrics
%=========environmental parameters======
cost.CO2NG=0.275; %carbon emission rate weight in kg/kw of NG
cost.CO2grid=0.469; %carbon emission rate weight in kg/kw of grid
cost.cCO2=1.2; %unit penalty cost of carbon emission
cost.CO2PV=0.0110;
cost.CO2Wind=0.02;
cost.CO2HP=0.100;
cost.CO2FC=0.2;
cost.CO2Elec=0.9;
cost.CO2bat=0.0111;

%% start of algorithm
scenario=1;% 1 2 3 negative trend positive
if scenario==1
     load('EV_negative_h.mat')
     EVCS_n=1;
elseif scenario==1
     load('EV_trend_h.mat')
     EVCS_n=2;
else scenario==2
    load('EV_Positive_h.mat')
    EVCS_n=5;
end

EVCS_max=2;
for y = 1:opt.Period
    Grid.nEV= size(EV_beh{y, 1}.v_is_charging,1);
    Charging_profile=double(EV_beh{y, 1}.v_is_charging);
    Grid.nEV_Private=floor(Grid.nEV);
    Grid.nEV_Public=ceil(Grid.nEV);
    Grid.EV_load_public= Charging_profile(1:Grid.nEV_Public,1:opt.Horizon);
    %Grid.EV_load_public=sum(Grid.EV_load_public,1);
    %Grid.EV_load_public(Grid.EV_load_public>22)=22;
    %[EVCS_state,waiting_time,Grid.EV_load]=getSimulationvalues1(Grid.nEV,EVCS_max,EVCS_n,Grid.EV_load_public,opt.Horizon) ;    
      [EVCS_state,waiting_time,Grid.EV_load]=getSimulationvalues(Grid.nEV,EVCS_max,EVCS_n,Grid.EV_load_public,opt.Horizon) ;
    busno=[2:ANZ.K];
    take_me = datasample(busno,Grid.nEV_Private);
    take_me=[4];
    Grid.EV_load_private= Charging_profile(:,1:opt.Horizon).*11;
    Grid.Pload(take_me,:)=Grid.Pload(take_me,:)+Grid.EV_load_private;
    Results = Deterministic_main(y,cost,Grid,Bat,heat,CHP,opt,ANZ,NET,mpc);
    %Results = Deterministic_main1(y,cost,Grid,Bat,heat,CHP,opt,ANZ,NET,mpc);
 end
End_time=toc


