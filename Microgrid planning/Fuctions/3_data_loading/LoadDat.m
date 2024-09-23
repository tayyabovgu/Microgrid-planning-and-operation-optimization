%% LoadDat.m

%% load grid data
% Nodes
[~,~,Gesamt]=xlsread(path.net_xlsx,'Nodes');                                                                       
NODE.StatL=Gesamt(2:end,1);                                     % station long name
NODE.StatB=Gesamt(2:end,2);                                     % station description
NODE.KNam=Gesamt(2:end,3);                                      % node name
NODE.KTyp=cell2mat(Gesamt(2:end,4));                            % node type
NODE.uKn=cell2mat(Gesamt(2:end,5));                             % nominal voltage in V
ANZ.K=size(NODE.KNam,1);                                        % number of nodes
% Lines
[~,~,Gesamt]=xlsread(path.net_xlsx,'Lines');                                                              
LIN.LtgNam=Gesamt(2:end,1);                                    % line name
LIN.vKNamHilf=Gesamt(2:end,2);                                 % name from node
LIN.zKNamHilf=Gesamt(2:end,3);                                 % name to node
LIN.vKStatHilf=Gesamt(2:end,4);                                % status from node (disconnector)
LIN.zKStatHilf=Gesamt(2:end,5);                                % status to node (disconnector)
LIN.StatLS=Gesamt(2:end,6);                                    % status of circuit breakers
LIN.StatE=Gesamt(2:end,7);                                     % status of earth connection
LIN.l=cell2mat(Gesamt(2:end,8));                               % line length in km
LIN.r=cell2mat(Gesamt(2:end,9));                               % resistance in Ohm/km
LIN.x=cell2mat(Gesamt(2:end,10));                              % reactance in Ohm/km
LIN.c=cell2mat(Gesamt(2:end,11));                              % capacitance in nF/km
LIN.Imax=cell2mat(Gesamt(2:end,13));                           % maximum current in kA
ANZ.L=size(LIN.LtgNam,1);                                      % number of lines                                              
% Transformers
[~,~,Gesamt]=xlsread(path.net_xlsx,'Trafos');                                                                      
TRA.TrNam=Gesamt(2:end,1);                                     % transformer name
TRA.vKNamHilf=Gesamt(2:end,2);                                 % name from node
TRA.zKNamHilf=Gesamt(2:end,3);                                 % name to node
TRA.vKStatHilf=Gesamt(2:end,4);                                % status from node (disconnector)
TRA.zKStatHilf=Gesamt(2:end,5);                                % status to node (disconnector)
TRA.StatLS=Gesamt(2:end,6);                                    % status of circuit breakers
TRA.StatE=Gesamt(2:end,7);                                     % status of earth connection
TRA.UrTOS=cell2mat(Gesamt(2:end,8));                           % primary voltage in kV
TRA.UrTUS=cell2mat(Gesamt(2:end,9));                           % secondary voltage in kV
TRA.SrT=cell2mat(Gesamt(2:end,10));                            % nominal apparent power in MVA
TRA.uk=cell2mat(Gesamt(2:end,11));                             % relative short circuit voltage in %
TRA.Pk=cell2mat(Gesamt(2:end,12));                             % short circuit losses in kW
TRA.P0=cell2mat(Gesamt(2:end,13));                             % no-load losses in kW
TRA.I0=cell2mat(Gesamt(2:end,14));                             % no-load current in %
ANZ.Tr=size(TRA.TrNam,1);                                      % number of transformers                                                                                
%Generation
[~,~,Gesamt]=xlsread(path.net_xlsx,'PV_generation');           % read generation data
GEN.GenNam=Gesamt(2:end,1);                                    % generation name
GEN.KNam=Gesamt(2:end,2);                                      % node name
GEN.PNenn=cell2mat(Gesamt(2:end,3));                           % active power in MW
GEN.Q=cell2mat(Gesamt(2:end,4));                               % reactive power in Mvar
GEN.ProfNam=Gesamt(2:end,5);                                   % Profile Name
ANZ.E=size(GEN.GenNam,1);                                      % number of generations                                                      
% Loads
[~,~,Gesamt]=xlsread(path.net_xlsx,'Loads');                                                                          
LOAD.LaNam=Gesamt(2:end,1);                                    % load name
LOAD.KNam=Gesamt(2:end,2);                                     % node name
LOAD.P=cell2mat(Gesamt(2:end,3));                              % active power in W
LOAD.Q=cell2mat(Gesamt(2:end,4));
ANZ.La=size(LOAD.LaNam,1);  
% LOAD.EV=cell2mat(Gesamt(2:end,12));
% ANZ.EV=sum(LOAD.EV);
 % Switches
        % Switches
        [~,~,Gesamt]=xlsread(path.net_xlsx,'Schalter');                     
        if size(Gesamt,1)>1                                                
            SWT.SwNam=Gesamt(2:end,1);                                     % switch name
            SWT.vKNamHilf=Gesamt(2:end,2);                                 % name from node
            SWT.zKNamHilf=Gesamt(2:end,3);                                 % name to node
            SWT.vKStatHilf=Gesamt(2:end,4);                                % status from node
            SWT.zKStatHilf=Gesamt(2:end,5);                                % status to node
            SWT.X=cell2mat(Gesamt(2:end,6));                               % reactance in Ohm
            SWT.Stat=cell2mat(Gesamt(2:end,7));                            % switching state
            ANZ.S=size(SWT.SwNam,1);                                       % number of switches
        else                                                                
            ANZ.S=0;                                                       % number of switches is equal to 0
            SWT=struct;
        end  
%%--------------------------------------------------------------------------------------------------------------------
%% load profile data
%%------------------------------------------------------------------------------------------------------------------------
% generation (P)
[~,~,Gesamt]=xlsread(path.prof_xlsx,'GEN_P');                    % active power generation
PROF.PGenProfNam=Gesamt(1,1:end);                                % generation names
PROF.PGenProf=Gesamt(2:end,1:end);                               % generation profiles
% generation (Q)
[~,~,Gesamt]=xlsread(path.prof_xlsx,'GEN_Q');                     % reactive power generation
PROF.QGenProfNam=Gesamt(1,1:end);                             % generation names
PROF.QGenProf=Gesamt(2:end,1:end);                            % generation profiles
% load (P)
[~,~,Gesamt]=xlsread(path.prof_xlsx,'LOAD_P');                % active power load
PROF.PLaNam=Gesamt(1,1:end);                                  % load names
PROF.PLaProfInd=Gesamt(2,1:end);% load profiles
% load (Q)
[~,~,Gesamt]=xlsread(path.prof_xlsx,'LOAD_Q');                % reactive power load
PROF.QLaNam=Gesamt(1,1:end);                                  % load names
PROF.QLaProfInd=Gesamt(2,1:end);                              % load profiles
%% load profiles based on HTW database
profile=load('LoadProfiles.mat');                             % load profile file
PROF.LProfsInit=profile.profile;% resave load profiles
clearvars Gesamt;                                                       % load not necessary variables
%% Allgemeine Parameter bestimmen
ANZ.tmax=size(PROF.PGenProf,1);                                     % number of time steps
        ANZ.Te=ANZ.L+ANZ.Tr+ANZ.S;                                          % number of terminals

    
