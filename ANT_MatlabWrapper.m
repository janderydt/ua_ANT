function ANT_MatlabWrapper(pgid,type)

UserVar.type = type;

logfile = pwd+"/jobs_master.log";
fid = fopen(logfile,'a+');

UserVar.fid = fid;

%% setup Ua folders
[~,hostname]= system("hostname"); 
if strfind(hostname,"C23000100")
    run /mnt/md0/Ua/setup_Ua2D.m
    UserVar.hostname = "C23000100";
    NumWorkers = 25;
end
UserVar.Table = pwd+"/RunTable.csv";

addpath(getenv("froot_tools"));

%% read table
RunTable = ANT_ReadWritetable(UserVar,[],'read');

Iexisting = find(RunTable{:,'ExpID'}~=0);
Inew = find(RunTable{:,'ExpID'}==0);

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
                "restart is required...",string(RunTable{ind,'ExpID'}));
            UserVar.Finished = 0;
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

                fprintf(fid,'============================\n');
                fprintf(fid,string(datetime("now"))+"\n");
                fprintf(fid,'============================\n');
                fprintf(fid,"> %s: Submitted.\n",UserVar.Experiment);
                
                something_submitted=1;
                Inew = [];

                % run info
                UserVar = ANT_GetUserVar_Diagnostic(RunTable,ind,UserVar);

                % launch Ua job
                UserVar = ANT_UaJob(RunTable,ind,UserVar,pgid);
    
            elseif type=="Inverse"
                
                while ~UserVar.Finished

                    something_submitted=1;
                    Inew = [];

                    % read Runtable again in case any changes were made by other
                    % processes
                    RunTable=ANT_ReadWritetable(UserVar,[],'read');
                    ind = find(RunTable{:,'ExpID'}(:) == UserVar.ExpID);
     
                    % initialize User variables
                    UserVar = ANT_GetUserVar_Inverse(RunTable,ind,UserVar);
    
                    % cumulative sum of number of iterations at the end of each
                    % inverse cycle
                    it_tmp = cumsum(UserVar.Inverse.Iterations);
    
                    %% Inverse cycle
                    if UserVar.InverseCycle
                    
                        while UserVar.Inverse.IterationsDone < it_tmp(UserVar.Inverse.Cycle) && ~UserVar.Finished
        
                            UserVar.TargetIterations = min(5000,it_tmp(UserVar.Inverse.Cycle)-UserVar.Inverse.IterationsDone);
    
                            UserVar = ANT_UaJob(RunTable,ind,UserVar,pgid);
    
                            %adjust Runtable
                            RunTable=ANT_ReadWritetable(UserVar,[],'read');
                            ind = find(RunTable{:,'ExpID'}(:) == UserVar.ExpID);
                            RunTable{ind,"InverseIterationsDone"} = UserVar.Inverse.IterationsDone;   
                            [~] = ANT_ReadWritetable(UserVar,RunTable,'write');
                
                        end   

                        fprintf(fid,'============================\n');
                        fprintf(fid,string(datetime("now"))+"\n");
                        fprintf(fid,'============================\n');
                        fprintf(UserVar.fid,"> %s: End inverse cycle %s.\n",UserVar.Experiment,string(UserVar.Inverse.Cycle));
    
                    %% Spinup cycle
                    elseif UserVar.SpinupCycle
    
                        UserVar = ANT_UaJob(RunTable,ind,UserVar,pgid);
    
                        %adjust Runtable
                        RunTable=ANT_ReadWritetable(UserVar,[],'read');
                        ind = find(RunTable{:,'ExpID'}(:) == UserVar.ExpID);
                        RunTable{ind,"SpinupYearsDone"} = UserVar.Spinup.YearsDone;   
                        [~] = ANT_ReadWritetable(UserVar,RunTable,'write');

                        fprintf(fid,'============================\n');
                        fprintf(fid,string(datetime("now"))+"\n");
                        fprintf(fid,'============================\n');
                        fprintf(UserVar.fid,"> %s: End spinup cycle %s.\n",UserVar.Experiment,string(UserVar.Spinup.Cycle));
    
                    end
        
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

    % generate unique ExpID and save to run table
    existingID = RunTable{:,"ExpID"};
    ExpID = 0;
    while ismember(ExpID,existingID) 
        ExpID = randi([1000 1999]);
    end
    UserVar.ExpID = ExpID;
    RunTable{ind,"ExpID"} = ExpID;
    RunTable{ind,'pgid'} = pgid;

    [~]=ANT_ReadWritetable(UserVar,RunTable,'write');

    % make copy of master folder for new experiment
    copyfile(['./ANT_',char(type),'_9999/'],['./ANT_',char(type),'_',num2str(ExpID),'/']); 

    % initialize some UserVars
    UserVar.Finished = 0;
    UserVar.Restart = 0;
    UserVar.Experiment = ['ANT_',char(type),'_',num2str(ExpID)];
    
    % now gather run info
    if type=="Diagnostic"

        fprintf(fid,'============================\n');
        fprintf(fid,string(datetime("now"))+"\n");
        fprintf(fid,'============================\n');
        fprintf(fid,"> %s: Submitted.\n",UserVar.Experiment);

        UserVar = ANT_GetUserVar_Diagnostic(RunTable,ind,UserVar);

        UserVar.Restart = 0;
    
        UserVar = ANT_UaJob(RunTable,ind,UserVar,pgid);

    elseif type=="Inverse"

        while ~UserVar.Finished
            
            % read Runtable again in case any changes were made by other
            % processes
            RunTable=ANT_ReadWritetable(UserVar,[],'read');
            ind = find(RunTable{:,'ExpID'}(:) == UserVar.ExpID);
               
            % initialize User variables
            UserVar = ANT_GetUserVar_Inverse(RunTable,ind,UserVar);
    
            % cumulative sum of number of iterations at the end of each
            % inverse cycle
            it_tmp = cumsum(UserVar.Inverse.Iterations);
    
            %% Inverse cycle
            if UserVar.InverseCycle
            
                while UserVar.Inverse.IterationsDone < it_tmp(UserVar.Inverse.Cycle) && ~UserVar.Finished
    
                    UserVar.TargetIterations = min(5000,it_tmp(UserVar.Inverse.Cycle)-UserVar.Inverse.IterationsDone);
    
                    UserVar = ANT_UaJob(RunTable,ind,UserVar,pgid);
    
                    %adjust Runtable
                    RunTable=ANT_ReadWritetable(UserVar,[],'read');
                    ind = find(RunTable{:,'ExpID'}(:) == UserVar.ExpID);
                    RunTable{ind,"InverseIterationsDone"} = UserVar.Inverse.IterationsDone;   
                    [~] = ANT_ReadWritetable(UserVar,RunTable,'write');

                end     

                fprintf(fid,'============================\n');
                fprintf(fid,string(datetime("now"))+"\n");
                fprintf(fid,'============================\n');
                fprintf(UserVar.fid,"> %s: End inverse cycle %s.\n",UserVar.Experiment,string(UserVar.Inverse.Cycle));
    
            %% Spinup cycle
            elseif UserVar.SpinupCycle
    
                UserVar = ANT_UaJob(RunTable,ind,UserVar,pgid);

                %adjust Runtable
                RunTable=ANT_ReadWritetable(UserVar,[],'read');
                ind = find(RunTable{:,'ExpID'}(:) == UserVar.ExpID);
                RunTable{ind,"SpinupYearsDone"} = UserVar.Spinup.YearsDone;   
                [~] = ANT_ReadWritetable(UserVar,RunTable,'write');

                fprintf(fid,'============================\n');
                fprintf(fid,string(datetime("now"))+"\n");
                fprintf(fid,'============================\n');
                fprintf(UserVar.fid,"> %s: End spinup cycle %s.\n",UserVar.Experiment,string(UserVar.Spinup.Cycle));
    
            end

        end
    end

    ANT_CleanUp(UserVar);

    diary off

else

    fprintf(fid,"   ...ANT_MatlabWrapper: Nothing to do. Try again later.\n");

end

