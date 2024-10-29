function ANT_CleanUp(UserVar)

%% This function only writes to experiment run tables

% read table
RunTable=ANT_ReadWritetable(UserVar,UserVar.runtable_exp,[],'read'); 
ind = find(RunTable{:,'ExpID'}==UserVar.ExpID); % all going well, this index should always be 1

% add flags and timestamps to table
if UserVar.Finished==1 && UserVar.Error==0

    RunTable{ind,"Submitted"} = 0;
    RunTable{ind,"Running"} = 0;
    RunTable{ind,"Finished"} = 1;
    RunTable{ind,"Restart"} = 0;
    RunTable{ind,"FinishedTime"} = string(datetime("now","format","yyyy-MM-dd HH:mm:ss"));
 
    fprintf(UserVar.fid_experimentlog,'> ANT_CleanUp: ExpID %s SUCCESSFULLY FINISHED.\n',string(UserVar.ExpID));    
    fprintf(UserVar.fid_experimentlog,'============================\n');

elseif UserVar.Finished==0 && UserVar.Restart==1

    RunTable{ind,"Submitted"} = 0;
    RunTable{ind,"Running"} = 0;
    RunTable{ind,"Finished"} = 0;
    RunTable{ind,"Restart"} = 1;
    RunTable{ind,"FinishedTime"} = string(datetime("now","format","yyyy-MM-dd HH:mm:ss"));
   
    fprintf(UserVar.fid_experimentlog,'> ANT_CleanUp: ExpID %s RESTART REQUIRED.\n',string(UserVar.ExpID));
    fprintf(UserVar.fid_experimentlog,'============================\n');

elseif UserVar.Finished==0 && UserVar.Restart==0 && UserVar.Error==0
	
    RunTable{ind,"Submitted"} = 0;
    RunTable{ind,"Running"} = 0;
    RunTable{ind,"Finished"} = 0;
    RunTable{ind,"Restart"} = 0;
    RunTable{ind,"FinishedTime"} = string(datetime("now","format","yyyy-MM-dd HH:mm:ss"));

    if isfield(UserVar,'InverseCycle')
        if UserVar.InverseCycle==1
	        fprintf(UserVar.fid_experimentlog,'> ANT_CleanUp: ExpID %s FINISHED INVERSE CYCLE.\n',string(UserVar.ExpID));
        else
            fprintf(UserVar.fid_experimentlog,'> ANT_CleanUp: ExpID %s FINISHED SPINUP CYCLE.\n',string(UserVar.ExpID));
        end
    end
    fprintf(UserVar.fid_experimentlog,'============================\n');

elseif UserVar.Error==1

    RunTable{ind,"Running"} = 0;
    RunTable{ind,"Error"} = 1;
    RunTable{ind,"ErrorTime"} = string(datetime("now","format","yyyy-MM-dd HH:mm:ss"));
   
    fprintf(UserVar.fid_experimentlog,'> ANT_CleanUp: ExpID %s ABORTED AND DID NOT FINISH.\n',string(UserVar.ExpID));    
    fprintf(UserVar.fid_experimentlog,'============================\n');

else

    error("ANT_Cleanup ExpID "+string(UsserVar.ExpID)+": unknown case");
    
end

RunTable{ind,"pgid"} = 0;

[~]=ANT_ReadWritetable(UserVar,UserVar.runtable_exp,RunTable,'write');
