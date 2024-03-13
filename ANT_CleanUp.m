function ANT_CleanUp(UserVar)

Table="RunTable.csv";

% read table
RunTable=ANT_ReadWritetable(UserVar,[],'read'); 
ind = find(RunTable{:,'ExpID'}==UserVar.ExpID);

% add flags and timestamps to table
if UserVar.Finished==1

    RunTable{ind,"Submitted"} = 0;
    RunTable{ind,"Running"} = 0;
    RunTable{ind,"Finished"} = 1;
    RunTable{ind,"Restart"} = 0;
    RunTable{ind,"FinishedTime"}(:) = datestr(now);

    fprintf(UserVar.fid,'============================\n');
    fprintf(UserVar.fid,string(datetime("now"))+"\n");
    fprintf(UserVar.fid,'============================\n');
    fprintf(UserVar.fid,'ExpID %s SUCCESSFULLY FINISHED.\n',string(UserVar.ExpID));    

else

    RunTable{ind,"Running"} = 0;
    RunTable{ind,"Error"} = 1;
    RunTable{ind,"ErrorTime"}(:) = datestr(now);
    
    fprintf(UserVar.fid,'============================\n');
    fprintf(UserVar.fid,string(datetime("now"))+"\n");
    fprintf(UserVar.fid,'============================\n');
    fprintf(UserVar.fid,'ExpID %s ABORTED AND DID NOT FINISH.\n',string(UserVar.ExpID));    
    
end

RunTable{ind,"pgid"} = 0;

[~]=ANT_ReadWritetable(UserVar,RunTable,'write');
