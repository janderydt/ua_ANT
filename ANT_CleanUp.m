function ANT_CleanUp(UserVar)

% read table
RunTable=ANT_ReadWritetable(UserVar,[],'read'); 
ind = find(RunTable{:,'ExpID'}==UserVar.ExpID);

% add flags and timestamps to table
if UserVar.Finished==1 && UserVar.Error==0

    RunTable{ind,"Submitted"} = 0;
    RunTable{ind,"Running"} = 0;
    RunTable{ind,"Finished"} = 1;
    RunTable{ind,"Restart"} = 0;
    RunTable{ind,"FinishedTime"}(:) = datestr(now);

    fprintf(UserVar.fid,'============================\n');
    fprintf(UserVar.fid,string(datetime("now"))+"\n");
    fprintf(UserVar.fid,'============================\n');
    fprintf(UserVar.fid,'ExpID %s SUCCESSFULLY FINISHED.\n',string(UserVar.ExpID));    

elseif UserVar.Finished==0 && UserVar.Restart==1

    RunTable{ind,"Submitted"} = 0;
    RunTable{ind,"Running"} = 0;
    RunTable{ind,"Restart"} = 1;

    fprintf(UserVar.fid,'============================\n');
    fprintf(UserVar.fid,string(datetime("now"))+"\n");
    fprintf(UserVar.fid,'============================\n');
    fprintf(UserVar.fid,'ExpID %s RESTART REQUIRED.\n',string(UserVar.ExpID)); 

elseif UserVar.Error==1

    RunTable{ind,"Running"} = 0;
    RunTable{ind,"Error"} = 1;
    RunTable{ind,"ErrorTime"}(:) = datestr(now);
    
    fprintf(UserVar.fid,'============================\n');
    fprintf(UserVar.fid,string(datetime("now"))+"\n");
    fprintf(UserVar.fid,'============================\n');
    fprintf(UserVar.fid,'ExpID %s ABORTED AND DID NOT FINISH.\n',string(UserVar.ExpID));    

else

    error("ANT_Cleanup: unknown case");
    
end

RunTable{ind,"pgid"} = 0;

[~]=ANT_ReadWritetable(UserVar,RunTable,'write');
