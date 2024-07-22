function UserVar = ANT_UaJob(RunTable,ind,UserVar,pgid)

%% copy latest Ua run files from master folder to experiment folder
copyfile("./ANT_"+UserVar.type+"_9999/*.m",UserVar.casefolder+"/"+UserVar.Experiment);

%% check if experiment runtable exists, otherwise, create a new one
if ~isfile(UserVar.runtable_exp)
    copyfile("./ANT_"+UserVar.type+"_9999/RunTable_ANT_"+UserVar.type+"_9999.csv",UserVar.runtable_exp);
end

%% log changes in tables
RunTable{ind,"SubmissionTime"}(:) = string(datetime("now")); 
RunTable{ind,"Submitted"} = 1;
RunTable{ind,"Running"} = 1;
RunTable{ind,'pgid'} = pgid;

% log changes in global table
%[~]=ANT_ReadWritetable(UserVar,UserVar.runtable_global,RunTable,'write');

% extract relevant row from global table
RunTable_exp = RunTable(ind,:);

% write to experiment runtable
[~]=ANT_ReadWritetable(UserVar,UserVar.runtable_exp,RunTable_exp,'write');

% launch job
cd(UserVar.casefolder+"/"+UserVar.Experiment);

% start diary
diary(UserVar.Experiment+".out")
diary on

UserVar.Error=0;

try
    
    UserVar = Ua2D(UserVar);   
    
catch ME

    fprintf(UserVar.fid_experimentlog,string(datetime("now"))+" || ERROR in the execution of ExpID %s.\n",string(UserVar.ExpID));

    UserVar.Breakout = 1;
    UserVar.Error = 1;

    msgString = getReport(ME,'extended');
    fprintf(UserVar.fid_experror,"%s \n\n",msgString);

end

if UserVar.Finished == 1
    UserVar.Breakout = 1;
end
        
cd ..
