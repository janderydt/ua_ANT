function ANT_GenerateRunTable_Diagnostic

% This script uses existing runtables for inverse simulations as a basis
% to construct a new (or add to an existing) runtable for diagnostic 
% simulations. Each inverse simulation is the starting point for 10 
% diagnostic simulations:
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
RunTable_inverse_files = "../ANT_Inverse/RunTable_ARCHER2_17-02-2025_"+string([24 25 26 27])+".csv";
Modeldomain = "AMUND";

%% year for original and new geometry
ExpID_oldgeom = 2000;
ExpID_newgeom = 2014; % 2009, 2014, 2018

%% which perturbations?
pert = ["Calv_thick"]; % Original, Calv, ISthick, thick, Calv_thick

%% construct new run table or append to existing run table
UserVar.Table = "./RunTable_ARCHER2_Diagnostic_"+Modeldomain+"_Umbi_"+string(ExpID_newgeom)+".csv";
if ~exist(UserVar.Table,"file")
    copyfile("EmptyTable.csv",UserVar.Table);
else
    copyfile(UserVar.Table,UserVar.Table+".orig");
end
UserVar.type = "Diagnostic";
RunTable_diagnostic = ANT_ReadWritetable(UserVar,UserVar.Table,[],'read');

for tt=1:numel(RunTable_inverse_files)
    %% read original run table
    UserVar.type = "Inverse";
    RunTable_inverse =  ANT_ReadWritetable(UserVar,RunTable_inverse_files(tt),[],'read');
          
    IndOriginalGeometry = find(contains(RunTable_diagnostic.Comments,"Original"));
    
    ExistingInverseRunsInDiagnosticTable=RunTable_diagnostic.InverseA(IndOriginalGeometry);

    if isempty(ExistingInverseRunsInDiagnosticTable)
        ExistingInverseRunsInDiagnosticTable = [];
    else
        if iscell(ExistingInverseRunsInDiagnosticTable)
            ExistingInverseRunsInDiagnosticTable = unique(cell2mat(ExistingInverseRunsInDiagnosticTable));
        else
            ExistingInverseRunsInDiagnosticTable = unique(ExistingInverseRunsInDiagnosticTable);
        end
    end
    
    BaseMesh=Modeldomain+"_"+string(ExpID_newgeom)+"_meshmin1500_meshmax100000_refined";
    
    for ii=1:size(RunTable_inverse,1)
        
        if ~ismember(RunTable_inverse.ExpID(ii),ExistingInverseRunsInDiagnosticTable) && RunTable_inverse.Finished(ii)
    
            Newrows=[];
            Newrow_tmp=[];
        
            for it=1:2
            
                Newrow_tmp = ...
                    {RunTable_inverse.Domain{ii},...    %Domain
                    0,...                               %pgid
                    0,...                               %ExpID
                    0,...                               %Submitted
                    "2000-01-01 00:00:00",...           %SubmissionTime
                    0,...                               %Running
                    0,...                               %Error
                    "2000-01-01 00:00:00",...           %ErrorTime
                    0,...                               %Finished
                    "2000-01-01 00:00:00",...           %FinishedTime
                    0,...                               %Restart
                    BaseMesh,...                        %BaseMesh
                    RunTable_inverse.ExpID(ii),...      %InverseA
        	        it,...                              %InverseCycleA
        	        RunTable_inverse.ExpID(ii),...      %InverseC
        	        it,...                              %InverseCycleC
        	        ExpID_oldgeom,...                   %Calv
        	        ExpID_oldgeom,...                   %ISthick
        	        it,...                              %InverseCycleIS (1=no spinup, 2=spinup)
        	        ExpID_oldgeom,...                   %GIthick
        	        it,...                              %InverseCycleGI (1=no spinup, 2=spinup)
                    ""                                  %Comments
                    };
        
                %% add details about geometry
                % Original geometry
                if ismember("Original",pert)
                    Newrow_tmp_orig = Newrow_tmp;
                    Newrow_tmp_orig{end} = "Original geometry "+string(RunTable_inverse.startGeometry(ii));
                    Newrow_tmp_orig{12} = strrep(BaseMesh,"_"+string(ExpID_newgeom)+"_","_"+string(ExpID_oldgeom)+"_");
                    RunTable_diagnostic = [RunTable_diagnostic ; Newrow_tmp_orig];
                end

                % Calv
                if ismember("Calv",pert)
                    Newrow_tmp_Calv = Newrow_tmp;
                    Newrow_tmp_Calv{17} = ExpID_newgeom;
                    Newrow_tmp_Calv{end} = "Ice front geometry "+string(ExpID_newgeom);
                    RunTable_diagnostic = [RunTable_diagnostic ; Newrow_tmp_Calv];
                end
        
                % ISthick
                if ismember("ISthick",pert)
                    Newrow_tmp_ISthick = Newrow_tmp;
                    Newrow_tmp_ISthick{12} = strrep(BaseMesh,"_"+string(ExpID_newgeom)+"_","_"+string(ExpID_oldgeom)+"_");
                    Newrow_tmp_ISthick{18} = ExpID_newgeom;
                    Newrow_tmp_ISthick{end} = "Ice shelf thickness "+string(ExpID_newgeom);
                    RunTable_diagnostic = [RunTable_diagnostic ; Newrow_tmp_ISthick];
                end
        
                % ISthick and GIthick      
                if ismember("thick",pert)
                    Newrow_tmp_thick = Newrow_tmp;
                    Newrow_tmp_thick{12} = strrep(BaseMesh,"_"+string(ExpID_newgeom)+"_","_"+string(ExpID_oldgeom)+"_");
                    Newrow_tmp_thick{18} = ExpID_newgeom;
                    Newrow_tmp_thick{20} = ExpID_newgeom;
                    Newrow_tmp_thick{end} = "Ice thickness "+string(ExpID_newgeom);
                    RunTable_diagnostic = [RunTable_diagnostic ; Newrow_tmp_thick];
                end
        
                % Calv, ISthick and GIthick
                if ismember("Calv_thick",pert)
                    Newrow_tmp_Calv_thick = Newrow_tmp;
                    Newrow_tmp_Calv_thick{17} = ExpID_newgeom; 
                    Newrow_tmp_Calv_thick{18} = ExpID_newgeom; 
                    Newrow_tmp_Calv_thick{20} = ExpID_newgeom; 
                    Newrow_tmp_Calv_thick{end} = "Ice front and thickness "+string(ExpID_newgeom);
                    RunTable_diagnostic = [RunTable_diagnostic ; Newrow_tmp_Calv_thick];
                end
    
            end
    
        end
    
    end
end

UserVar.type = "Diagnostic";
ANT_ReadWritetable(UserVar,UserVar.Table,RunTable_diagnostic,'write');