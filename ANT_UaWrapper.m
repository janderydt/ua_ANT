function ANT_UaWrapper(pgid,type)

if nargin==2 
    % ensure correct format
    if ischar(pgid)
        pgid = str2double(pgid);
    end
    if ischar(type)
        type = string(type);
    end
    walltime = 1e10; % set to some large number
else
    pgid=[];
    type=[];
    walltime=[];
end

%% initialize log file
UserVar.home = pwd;

logfile = UserVar.home+"/jobs_master.log";
fid = fopen(logfile,'a+');

UserVar.fid = fid;

%% find host and setup matlab path
[~,hostname]= system("hostname"); 
if strfind(hostname,"C23000100")
    run /mnt/md0/Ua/setup_Ua2D.m
    UserVar.hostname = "C23000100";
    NumWorkers = 25;
elseif strfind(hostname,"sauron")
    run /home/wchm8/Documents/Ua/setup_Ua2D.m
    UserVar.hostname = "sauron";
elseif strfind(hostname,"nid") % ARCHER2
    UserVar.hostname = "ARCHER2";
else
    error("Hostname "+hostname+" not found.");
end

if ~contains(UserVar.hostname,"ARCHER2")
    addpath(getenv("froot_tools"));
end

%% obtain run type and other inputs if not already provided
if nargin<2 
    if contains(UserVar.hostname,"ARCHER2")
        % on ARCHER2 there are no inputs. instead we read a text file with config variables
        configfile = pwd+"/ua_config.txt";
        if ~exist(configfile,"file")
            error("Specify config file"+configfile+".");
        end
        fileID = fopen(configfile,"r"); 
        tline = fgetl(fileID);
        while ischar(tline)	
            if ~isempty(tline)
                if ~startsWith(tline(1), '#')
                    % process line
                    if contains(tline,'runtype')
                        type = string(erase(tline,["runtype"," ","=",""""]));
                    elseif contains(tline,'pgid')
                        pgid = str2num(erase(tline,["pgid"," ","="]));
                    elseif contains(tline,'walltime')
                        walltime = seconds(duration(erase(tline,["walltime"," ","="])));
                    end
                end
            end
            % Now read the next line.
            tline = fgetl(fileID);
        end
        % All done reading all lines, so close the file.
        fclose(fileID);      
    else
        error("Running job on "+UserVar.hostname+" so need 2 inputs, got "+string(nargin)+" instead: pgid ("+string(pgid)+") and runtype ("+string(type)+").");
    end
end

%% check that all inputs are now available
if isempty(type) || isempty(walltime) || isempty(pgid)
    error("One of the following mandatory input variables is empty: type ("+string(type)+"), "+ ...
    "walltime ("+string(walltime)+"), pgid ("+string(pgid)+").");
else
    UserVar.type = type;
    UserVar.walltime = walltime-15*60; % subtract 10min for delays at the start; this is very common on ARCHER2
    UserVar.pgid = pgid;
end

%% read run table
UserVar.Table = pwd+"/RunTable_"+UserVar.hostname+".csv";

RunTable = ANT_ReadWritetable(UserVar,[],'read');

if ~isempty(RunTable)
    Iexisting = find(RunTable{:,'ExpID'}~=0);
    Inew = find(RunTable{:,'ExpID'}==0);
else
    fprintf(fid,"Empty RunTable - stop job./n");
    return
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
        UserVar.Domain = RunTable{ind,'Domain'};
        UserVar.Experiment = [char(UserVar.Domain),'_',char(type),'_',num2str(RunTable{ind,'ExpID'})];
        UserVar.ExpID = RunTable{ind,'ExpID'};
    
        % check if submitted but not running
        indsnr = RunTable{ind,'Submitted'}~=0 & RunTable{ind,'Running'}==0;
        % check if not submitted, not running, not finished
        indnsnr = RunTable{ind,'Submitted'}==0 & RunTable{ind,'Running'}==0 & RunTable{ind,'Finished'}==0;
    
        % something wrong?
        if indsnr
            fprintf(fid,"   ...ANT_UaWrapper: ExpID %s has been submitted, but corresponding jobID has not " + ...
                "been found. Either something went wrong, or the run has " + ...
                "finished. Check log files for errors.\n",string(RunTable{indsnr,'ExpID'}));
            error('');
        end

        
        if indnsnr
            fprintf(fid,"   ...ANT_UaWrapper: ExpID %s has been not yet been submitted. Let's check if a " + ...
                "restart is required...",string(RunTable{ind,'ExpID'}));

            % initialze some variables
            UserVar.Finished = 0;
            UserVar.Breakout = 0;

            % new run or restart?
            if RunTable{ind,'Restart'}==1
                UserVar.Restart = 1;
                fprintf(fid,"yes.\n");
                fprintf(fid,"   ...ANT_UaWrapper: Restarting ExpID %s...\n",string(RunTable{ind,'ExpID'}));
            else
                UserVar.Restart = 0;
                fprintf(fid,"no.\n");
            end

            % now gather run info and launch job
            if type=="Diagnostic"

                fprintf(fid,'============================\n');
                fprintf(fid,string(datetime("now"))+"\n");
                fprintf(fid,'============================\n');
                fprintf(fid,"> %s: Submitted on %s.\n",UserVar.Experiment,UserVar.hostname);
                
                something_submitted=1;
                Inew = [];

                % run info
                UserVar = ANT_GetUserVar_Diagnostic(RunTable,ind,UserVar);

                % launch Ua job
                UserVar = ANT_UaJob(RunTable,ind,UserVar,pgid);
    
            elseif type=="Inverse"

                while ~UserVar.Breakout

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
                    
                        while (UserVar.Inverse.IterationsDone < it_tmp(UserVar.Inverse.Cycle) && ~UserVar.Breakout)

                            if contains(UserVar.hostname,"ARCHER2")
                                nit = 10000;
                            else
                                nit = 5000;
                            end
        
                            UserVar.TargetIterations = min(nit,it_tmp(UserVar.Inverse.Cycle)-UserVar.Inverse.IterationsDone);
                            fprintf(UserVar.fid,"> Doing %s iterations.\n",string(UserVar.TargetIterations));
    
                            UserVar = ANT_UaJob(RunTable,ind,UserVar,pgid);
    
                            %adjust Runtable
                            RunTable=ANT_ReadWritetable(UserVar,[],'read');
                            ind = find(RunTable{:,'ExpID'}(:) == UserVar.ExpID);
                            RunTable{ind,"InverseIterationsDone"} = UserVar.Inverse.IterationsDone;   
                            RunTable{ind,"Restart"} = UserVar.Restart;
                            [~] = ANT_ReadWritetable(UserVar,RunTable,'write');
                
                        end   

                        it_tmp = cumsum(UserVar.Inverse.Iterations);

                        fprintf(fid,'============================\n');
                        fprintf(fid,string(datetime("now"))+"\n");
                        fprintf(fid,'============================\n');
                        fprintf(UserVar.fid,"> %s: Breaking out of inverse cycle %s. Done %s iterations out of %s.\n",...
                        UserVar.Experiment,string(UserVar.Inverse.Cycle),string(UserVar.Inverse.IterationsDone),string(it_tmp(end)));
    
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
        ExpID = randi([2000 2999]);
    end
    UserVar.ExpID = ExpID;
    RunTable{ind,"ExpID"} = ExpID;
    RunTable{ind,'pgid'} = pgid;

    [~]=ANT_ReadWritetable(UserVar,RunTable,'write');
    
    % initialize some UserVars
    UserVar.Finished = 0;
    UserVar.Restart = 0;
    UserVar.Breakout = 0;
    UserVar.Domain = RunTable{ind,'Domain'};
    UserVar.Experiment = [char(UserVar.Domain),'_',char(type),'_',num2str(ExpID)];

    % make copy of master folder for new experiment
    % if new folder already exists: rename first
    sourcefolder = ['./ANT_',char(type),'_9999/'];
    newfolder = ['./',char(UserVar.Domain),'_',char(type),'_',num2str(ExpID)];
    if exist(newfolder,"dir") == 7
        movefile(newfolder,[newfolder,'_old/']);
    else
        copyfile(sourcefolder,[newfolder,'/']); 
    end
    
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

        while ~UserVar.Breakout
            
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
            
                while (UserVar.Inverse.IterationsDone < it_tmp(UserVar.Inverse.Cycle) && ~UserVar.Breakout)
    
                    if contains(UserVar.hostname,"ARCHER2")
                        nit = 10000;
                    else
                        nit = 5000;
                    end

                    UserVar.TargetIterations = min(nit,it_tmp(UserVar.Inverse.Cycle)-UserVar.Inverse.IterationsDone);
                    fprintf(UserVar.fid,"> Doing %s iterations.\n",string(UserVar.TargetIterations));
    
                    UserVar = ANT_UaJob(RunTable,ind,UserVar,pgid);
    
                    %adjust Runtable
                    RunTable=ANT_ReadWritetable(UserVar,[],'read');
                    ind = find(RunTable{:,'ExpID'}(:) == UserVar.ExpID);
                    RunTable{ind,"InverseIterationsDone"} = UserVar.Inverse.IterationsDone;  
                    RunTable{ind,"Restart"} = UserVar.Restart;
                    [~] = ANT_ReadWritetable(UserVar,RunTable,'write');

                end     

                it_tmp = cumsum(UserVar.Inverse.Iterations);

                fprintf(fid,'============================\n');
                fprintf(fid,string(datetime("now"))+"\n");
                fprintf(fid,'============================\n');
                fprintf(UserVar.fid,"> %s: Breaking out of inverse cycle %s. Done %s iterations out of %s.\n",...
                    UserVar.Experiment,string(UserVar.Inverse.Cycle),string(UserVar.Inverse.IterationsDone),string(it_tmp(end)));
    
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
                fprintf(UserVar.fid,"> %s: Breaking out of spinup cycle %s.\n",UserVar.Experiment,string(UserVar.Spinup.Cycle));
    
            end

        end
    end

    ANT_CleanUp(UserVar);

    diary off

else

    fprintf(fid,"   ...ANT_UaWrapper: Nothing to do. Try again later.\n");

end
