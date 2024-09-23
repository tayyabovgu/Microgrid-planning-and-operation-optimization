function path = Allp(SetC)
%% Define Filepaths for Various Data Types
% SetC: Structure that holds configuration options, including the network case.

% Base directory (where functions are located)
path.dat1 = fullfile(pwd, 'Functions');
% Subdirectories for different data types
path.Grid_data              = fullfile(path.dat1, '1_Grid_files');               % Grid data
path.Load_profiles          = fullfile(path.dat1, '2_Load_Profiles');              % load profile
path.Load_dataloading       = fullfile(path.dat1, '3_data_loading');              % load profile
path.Load_interpolation     = fullfile(path.dat1, '4_Interpolation');              % load profile
path.Load_interpolation     = fullfile(path.dat1, '5_Loadflow');              % load profile
path.Load_interpolation     = fullfile(path.dat1, '6_Tarrif');              % load profile
path.Load_heating     = fullfile(path.dat1, '7_Heating');              % load profile
% % Specific files that depend on the selected network case (SetC.netcase)
path.net_xlsx=fullfile(path.Grid_data,[SetC.netcase,'.xlsx']);
path.prof_xlsx=fullfile(path.Load_profiles,[SetC.netcase,'_prof.xlsx']);
path.heat_xlsx=fullfile(path.Load_heating,[SetC.netcase,'_heatgrid.xlsx']);
end
