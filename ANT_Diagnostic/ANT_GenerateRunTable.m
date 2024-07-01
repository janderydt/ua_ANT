function ANT_GenerateRunTable

% This script uses an existing runtable for inverse simulations as a basis
% to construct a new runtable for diagnostic simulations. Each inverse
% simulation is the starting point for 10 diagnostic simulations:
% 1. a uv solve with the original configuration (with/without spin-up)
% 2. a uv solve with a new calving front (with/without spin-up)
% 3. a uv solve with a new ice shelf thickness (with/without spin-up)
% 4. a uv solve with a new ice shelf and grounded ice thickness 
% with/without spin-up)
% 5. a uv solve with a new calving front, ice shelf thickness and grounded
% ice thickness (with/without spin-up)
% NOTE: The A and C fields remain unchanged in these experiments.

addpath("../");
UserVar.home = pwd;

%% original inversion run table
RunTable_inverse_file = "../ANT_Inverse/RunTable_ARCHER2_2.csv";

%% experiment number for new geometry
ExpID_newgeom = 1111;

%% read original run table
UserVar.type = "Inverse";
RunTable_inverse =  ANT_ReadWritetable(UserVar,RunTable_inverse_file,[],'read');

%% construct new run table
UserVar.Table = erase(RunTable_inverse_file,["../ANT_Inverse/",".csv"])+"_Diagnostic.csv";
if ~exist(UserVar.Table,"file")
    copyfile("EmptyTable.csv",UserVar.Table);
end
UserVar.type = "Diagnostic";
RunTable_diagnostic = ANT_ReadWritetable(UserVar,UserVar.Table,[],'read');
BaseMesh='2000_2009_2014_2018_meshmin3000_meshmax100000_refined';

for ii=1:size(RunTable_inverse,1)
    
    Newrows=[];
    Newrow_tmp=[];

    for it=1:2
    
        Newrow_tmp = ...
            {RunTable_inverse.Domain{ii},...  %Domain
            0,...                               %pgid
            0,...                               %ExpID
            0,...                               %Submitted
            "01/01/2000 00:00:00",...           %SubmissionTime
            0,...                               %Running
            0,...                               %Error
            "01/01/2000 00:00:00",...           %ErrorTime
            0,...                               %Finished
            "01/01/2000 00:00:00",...           %FinishedTime
            0,...                               %Restart
            BaseMesh,...                        %BaseMesh
            RunTable_inverse.ExpID(ii),...      %InverseA
        	RunTable_inverse.ExpID(ii),...      %InverseAFill (fields will be extruded)
        	it,...                              %InverseCycleA
        	RunTable_inverse.ExpID(ii),...      %InverseC
        	RunTable_inverse.ExpID(ii),...      %InverseCFill (fields will be extruded)
        	it,...                              %InverseCycleC
        	RunTable_inverse.ExpID(ii),...      %Calv
        	RunTable_inverse.ExpID(ii),...      %ISthick
        	it,...                              %InverseCycleIS (1=no spinup, 2=spinup)
        	RunTable_inverse.ExpID(ii),...      %GIthick
        	it,...                              %InverseCycleGI (1=no spinup, 2=spinup)
            ""                                  %Comments
            };

        %% add details about geometry
        % none
        Newrow_tmp_orig = Newrow_tmp;
        Newrow_tmp_orig{end} = "Original geometry "+string(RunTable_inverse.startGeometry(ii));
        RunTable_diagnostic = [RunTable_diagnostic ; Newrow_tmp_orig];

        % Calv
        Newrow_tmp_Calv = Newrow_tmp;
        Newrow_tmp_Calv{19} = ExpID_newgeom;
        Newrow_tmp_Calv{end} = "Ice front geometry 2018";
        RunTable_diagnostic = [RunTable_diagnostic ; Newrow_tmp_Calv];

        % ISthick
        Newrow_tmp_ISthick = Newrow_tmp;
        Newrow_tmp_ISthick{20} = ExpID_newgeom;
        Newrow_tmp_ISthick{end} = "Ice shelf thickness 2018";
        RunTable_diagnostic = [RunTable_diagnostic ; Newrow_tmp_ISthick];

        % ISthick and GIthick       
        Newrow_tmp_ISthick_GIthick = Newrow_tmp;
        Newrow_tmp_ISthick_GIthick{20} = ExpID_newgeom;
        Newrow_tmp_ISthick_GIthick{22} = ExpID_newgeom;
        Newrow_tmp_ISthick_GIthick{end} = "Ice thickness 2018";
        RunTable_diagnostic = [RunTable_diagnostic ; Newrow_tmp_ISthick_GIthick];

        % Calv, ISthick and GIthick
        Newrow_tmp_Calv_ISthick_GIthick = Newrow_tmp;
        Newrow_tmp_Calv_ISthick_GIthick{19} = ExpID_newgeom; 
        Newrow_tmp_Calv_ISthick_GIthick{20} = ExpID_newgeom; 
        Newrow_tmp_Calv_ISthick_GIthick{22} = ExpID_newgeom; 
        Newrow_tmp_Calv_ISthick_GIthick{end} = "Ice front and thickness 2018";
        RunTable_diagnostic = [RunTable_diagnostic ; Newrow_tmp_Calv_ISthick_GIthick];
    end

end

ANT_ReadWritetable(UserVar,UserVar.Table,RunTable_diagnostic,'write');