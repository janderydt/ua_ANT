function ANT_UaJob(RunTable,ind,UserVar,pgid,fid)

% log changes in table
RunTable{ind,"SubmissionTime"}(:) = datestr(now); 
RunTable{ind,"Submitted"} = 1;
RunTable{ind,"Running"} = 1;
RunTable{ind,'pgid'} = pgid;
writetable(RunTable,Table);

% launch job
cd(UserVar.Experiment);

diary(UserVar.Experiment+".out")
diary on

something_submitted = 1;

fprintf(fid,'============================\n');
fprintf(fid,string(datetime("now"))+"\n");
fprintf(fid,'============================\n');
fprintf(fid,"> ANT_MatlabWrapper: (Re-)Submitted %s.\n",UserVar.Experiment);

Inew = [];

try
    UserVar = Ua2D(UserVar);                
catch ME
    fprintf(fid,'============================\n');
    fprintf(fid,string(datetime("now"))+"\n");
    fprintf(fid,'============================\n');
    fprintf(fid,'An error occurred in the execution of ExpID %s.\n',string(UserVar.ExpID));

    UserVar.Finished = 0;

    msgString = getReport(ME,'extended'); 
    fprintf(fid,"%s \n",msgString);
end

if UserVar.Finished
    UserVar.Restart = 1;
end
        
cd ..

% read Runtable again in case any changes were made by other
% processes
RunTable=readtable(Table); 
ind = find(RunTable{:,'ExpID'}(:) == UserVar.ExpID);     
writetable(RunTable,Table);