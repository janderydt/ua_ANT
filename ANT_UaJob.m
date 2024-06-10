function UserVar = ANT_UaJob(RunTable,ind,UserVar,pgid)

%% copy latest Ua run files from master folder to experiment folder
copyfile("./ANT_"+UserVar.type+"_9999/*.m",UserVar.Experiment);

%% check if experiment runtable exists, otherwise, create a new one
if ~isfile(UserVar.runtable_exp)
    copyfile("./ANT_Inverse_9999/RunTable_ANT_Inverse_9999.csv",UserVar.runtable_exp);
end

%% log changes in experiment table
% modify content of global table
RunTable{ind,"SubmissionTime"}(:) = string(datetime("now")); 
RunTable{ind,"Submitted"} = 1;
RunTable{ind,"Running"} = 1;
RunTable{ind,'pgid'} = pgid;

% extract relevant row from global table
RunTable_exp = RunTable(ind,:);

% write to experiment runtable
[~]=ANT_ReadWritetable(UserVar,UserVar.runtable_exp,RunTable_exp,'write');

% launch job
cd(UserVar.Experiment);

% start diary
diary(UserVar.Experiment+".out")
diary on

something_submitted = 1;

Inew = [];
UserVar.Error=0;

try
    fid = fopen( UserVar.home+"/ua_submitted", 'wt' );
    fclose(fid);

    UserVar = Ua2D(UserVar);   
    
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

if UserVar.Finished == 1
    UserVar.Breakout = 1;
end
        
cd ..
