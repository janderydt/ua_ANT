function ANT_MatlabWrapper(pgid,type)

cwd = pwd;

logfile = cwd+"/jobs_master.log";
fid = fopen(logfile,'a+');

UserVar.fid = fid;

%% setup Ua folders
[~,hostname]= system("hostname"); 
if strfind(hostname,"C23000100")
    run /mnt/md0/Ua/setup_Ua2D.m
    UserVar.hostname = "C23000100";
    NumWorkers = 25;
end
Table = cwd+"/RunTable.csv";

addpath(getenv("froot_tools"));

%% read table
if exist(Table,'file')
    RunTable=readtable(Table); 
    Iexisting = find(RunTable{:,'ExpID'}~=0);
    Inew = find(RunTable{:,'ExpID'}==0);
else    
    error("Runtable does not exist"); 
end

%% launch jobs
UserVar.Restart = 0;
something_submitted = 0; kk=1;

% first deal with existing jobs
if ~isempty(Iexisting)

    while something_submitted==0 && kk<=numel(Iexisting)
   
        ind = Iexisting(kk);

        UserVar.Experiment = ['ANT_',char(type),'_',num2str(RunTable{ind,'ExpID'})];
        UserVar.ExpID = RunTable{ind,'ExpID'};
    
        % check submitted but not running
        indsnr = RunTable{ind,'Submitted'}~=0 & RunTable{ind,'Running'}==0;
        % check not submitted, not running, not finished
        indnsnr = RunTable{ind,'Submitted'}==0 & RunTable{ind,'Running'}==0 & RunTable{ind,'Finished'}==0;
    
        % something wrong?
        if indsnr
            fprintf(fid,"   ...ANT_MatlabWrapper: ExpID %s has been submitted, but corresponding jobID has not " + ...
                "been found. Either something went wrong, or the run has" + ...
                "finished. Let's check.../n",string(RunTable{indsnr,'ExpID'}));
            % TO DO: check if run has finished

        elseif indnsnr
            fprintf(fid,"   ...ANT_MatlabWrapper: ExpID %s has been not yet been submitted. Let's check if a " + ...
                "restart is required...\n",string(RunTable{ind,'ExpID'}));
            % TO DO: check run has finished
            if RunTable{ind,'Restart'}==1
                UserVar.Restart = 1;
                fprintf(fid,"   ...ANT_MatlabWrapper: Restarting ExpID %s...\n",string(RunTable{ind,'ExpID'}));
		        RunTable{ind,"Restart"} = 0;
            end

            % now gather run info
            if type=="Diagnostic"
                UserVar = ANT_GetUserVar_Diagnostic(RunTable,ind,UserVar);
            elseif type=="Inverse"
                UserVar = ANT_GetUserVar_Inverse(RunTable,ind,UserVar);
            end
            
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
            fprintf(fid,"> ANT_MatlabWrapper: (Re-)Submitted %s. \n",UserVar.Experiment);

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

            cd ..
    
            % read Runtable again in case any changes were made by other
            % processes
            RunTable=readtable(Table); 
            ind = find(RunTable{:,'ExpID'}(:) == UserVar.ExpID);
            writetable(RunTable,Table);

            ANT_CleanUp(UserVar);

            diary off
    
        end

        kk=kk+1;

    end
    
end

% deal with new jobs if all existing jobs have been dealt with
if ~isempty(Inew)

    ind = Inew(1);

    % generate unique ExpID and copy default scripts to new folder
    existingID = RunTable{:,"ExpID"};
    ExpID = 0;
    while ismember(ExpID,existingID) 
        ExpID = randi([1000 1999]);
    end
    UserVar.ExpID = ExpID;
    RunTable{ind,"ExpID"} = ExpID;
    RunTable{ind,'pgid'} = pgid;
    
    copyfile(['./ANT_',char(type),'_9999/'],['./ANT_',char(type),'_',num2str(ExpID),'/']); 

    % initialize UserVar.Finished
    UserVar.Finished = 1;

    % now gather run info
    UserVar.Experiment = ['ANT_',char(type),'_',num2str(ExpID)];
    if type=="Diagnostic"
        UserVar = ANT_GetUserVar_Diagnostic(RunTable,ind,UserVar);
    elseif type=="Inverse"
        UserVar = ANT_GetUserVar_Inverse(RunTable,ind,UserVar);
    end
    UserVar.Restart = 0;

    % adjust table
    RunTable{ind,"Submitted"} = 1;
    RunTable{ind,"Running"} = 1;
    RunTable{ind,"SubmissionTime"}(:) = datestr(now);        
    writetable(RunTable,Table);

    % launch job
    cd(UserVar.Experiment);

    diary(UserVar.Experiment+".out")
    diary on

    fprintf(fid,'============================\n');
    fprintf(fid,string(datetime("now"))+"\n");
    fprintf(fid,'============================\n');
    fprintf(fid,"> ANT_MatlabWrapper: (Re-)Submitted %s. \n",UserVar.Experiment);

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

    cd ..

    % read Runtable again in case any changes were made by other
    % processes
    RunTable=readtable(Table); 
    ind = find(RunTable{:,'ExpID'}(:) == UserVar.ExpID);
    writetable(RunTable,Table);   

    ANT_CleanUp(UserVar);

    diary off

else

    fprintf(fid,"   ...ANT_MatlabWrapper: Nothing to do. Try again later.\n");

end

