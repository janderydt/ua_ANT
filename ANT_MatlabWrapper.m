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
% initialize some variables
UserVar.Restart = 0;
something_submitted = 0; kk=1;

% first deal with existing jobs
if ~isempty(Iexisting)

    while something_submitted==0 && kk<=numel(Iexisting)
   
        % table row
        ind = Iexisting(kk);

        % experiment name and unique ID
        UserVar.Experiment = ['ANT_',char(type),'_',num2str(RunTable{ind,'ExpID'})];
        UserVar.ExpID = RunTable{ind,'ExpID'};
    
        % check if submitted but not running
        indsnr = RunTable{ind,'Submitted'}~=0 & RunTable{ind,'Running'}==0;
        % check if not submitted, not running, not finished
        indnsnr = RunTable{ind,'Submitted'}==0 & RunTable{ind,'Running'}==0 & RunTable{ind,'Finished'}==0;
    
        % something wrong?
        if indsnr
            fprintf(fid,"   ...ANT_MatlabWrapper: ExpID %s has been submitted, but corresponding jobID has not " + ...
                "been found. Either something went wrong, or the run has" + ...
                "finished. Check log files for errors./n",string(RunTable{indsnr,'ExpID'}));
            error('');
        end

        
        if indnsnr
            fprintf(fid,"   ...ANT_MatlabWrapper: ExpID %s has been not yet been submitted. Let's check if a " + ...
                "restart is required...\n",string(RunTable{ind,'ExpID'}));
            % new run or restart?
            if RunTable{ind,'Restart'}==1
                UserVar.Restart = 1;
                fprintf(fid,"yes.\n");
                fprintf(fid,"   ...ANT_MatlabWrapper: Restarting ExpID %s...\n",string(RunTable{ind,'ExpID'}));
            else
                UserVar.Restart = 0;
                fprintf(fid,"no.\n");
            end

            % now gather run info and launch job
            if type=="Diagnostic"

                % run info
                UserVar = ANT_GetUserVar_Diagnostic(RunTable,ind,UserVar);

                % launch Ua job
                ANT_UaJob(RunTable,ind,UserVar,pgid,fid)
    
            elseif type=="Inverse"
                
                % initialize variables
                UserVar = ANT_GetUserVar_Inverse(RunTable,ind,UserVar,fid);

                % cumulative sum of number of iterations at the end of each
                % inverse cycle
                it_tmp = cumsum(UserVar.Inverse.TargetIterations);

                %% Inverse cycle
                if UserVar.InverseCycle
                
                    while UserVar.Inverse.IterationsDone < it_tmp(UserVar.Inverse.Cycle) && UserVar.Finished
    
                        UserVar.Iterations = min(5000,it_tmp(UserVar.Inverse.Cycle)-UserVar.Inverse.IterationsDone);

                        ANT_UaJob(RunTable,ind,UserVar,pgid,fid);

                        RunTable{ind,"IterationsDone"} = UserVar.Inverse.IterationsDone;   
                        writetable(RunTable,Table);
            
                    end            

                %% Spinup cycle
                elseif UserVar.SpinupCycle

                    ANT_UaJob(RunTable,ind,UserVar,pgid,fid);

                end

            end

            % 
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

    % initialize UserVars
    UserVar.Finished = 1;
    UserVar.Restart = 0;

    % now gather run info
    UserVar.Experiment = ['ANT_',char(type),'_',num2str(ExpID)];

    if type=="Diagnostic"

        UserVar = ANT_GetUserVar_Diagnostic(RunTable,ind,UserVar);

        UserVar.Restart = 0;
    
        ANT_UaJob(RunTable,ind,UserVar,pgid,fid);


    elseif type=="Inverse"

        % cumulative sum of number of iterations at the end of each
        % inverse cycle
        it_tmp = cumsum(UserVar.Inverse.TargetIterations);

        %% Inverse cycle
        if UserVar.InverseCycle
        
            while UserVar.Inverse.IterationsDone < it_tmp(UserVar.Inverse.Cycle) && UserVar.Finished

                UserVar.Iterations = min(5000,it_tmp(UserVar.Inverse.Cycle)-UserVar.Inverse.IterationsDone);

                ANT_UaJob(RunTable,ind,UserVar,pgid,fid);

                RunTable{ind,"IterationsDone"} = UserVar.Inverse.IterationsDone;   
                writetable(RunTable,Table);
    
            end
    

        %% Spinup cycle
        elseif UserVar.SpinupCycle

            ANT_UaJob(RunTable,ind,UserVar,pgid,fid);

        end

    end

    ANT_CleanUp(UserVar);

    diary off

else

    fprintf(fid,"   ...ANT_MatlabWrapper: Nothing to do. Try again later.\n");

end

