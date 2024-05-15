function UserVar = ANT_UaJob(RunTable,ind,UserVar,pgid)

% log changes in table
RunTable{ind,"SubmissionTime"}(:) = datestr(now); 
RunTable{ind,"Submitted"} = 1;
RunTable{ind,"Running"} = 1;
RunTable{ind,'pgid'} = pgid;
[~]=ANT_ReadWritetable(UserVar,RunTable,'write');

% copy latest Ua run files from master folder to experiment folder
copyfile("./ANT_"+UserVar.type+"_9999/*.m",UserVar.Experiment);

% launch job
cd(UserVar.Experiment);

diary(UserVar.Experiment+".out")
diary on

something_submitted = 1;

Inew = [];

try
    UserVar = Ua2D(UserVar);   
    UserVar.Error = 0;
    
catch ME
    fprintf(UserVar.fid,'============================\n');
    fprintf(UserVar.fid,string(datetime("now"))+"\n");
    fprintf(UserVar.fid,'============================\n');
    fprintf(UserVar.fid,'An error occurred in the execution of ExpID %s.\n',string(UserVar.ExpID));

    UserVar.Breakout = 1;
    UserVar.Error = 1;

    msgString = getReport(ME,'extended'); 
    fprintf(UserVar.fid,"%s \n",msgString);
end

if UserVar.hostname == "ARCHER2"
    UserVar.Breakout = 1;
elseif UserVar.Finished == 1
    UserVar.Breakout = 1;
end
        
cd ..