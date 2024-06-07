function ANT_UaWrapper(ua_config,pgid,type)

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

%% initialize log file
UserVar.home = pwd;

logfile = UserVar.home+"/jobs_master_"+UserVar.hostname+".log";
fid = fopen(logfile,'a+');

UserVar.fid = fid;

%% deal with inputs

if nargin==3
    % ensure correct format
    ua_config = string(ua_config);
    if ischar(pgid)
        pgid = str2double(pgid);
    end
    if ischar(type)
        type = string(type);
    end
elseif nargin==1
    ua_config = string(ua_config);
    pgid=[];
    type=[];
else
    error("You need to either provide a config file or pgid and type. Instead got none.");
end
    
if ~isempty(ua_config)
    %% read inputs from config file
    configfile = UserVar.home+"/"+ua_config;
    if ~exist(configfile,"file")
        error("Config file "+configfile+" does not exist.");
    end
    fileID = fopen(configfile,"r"); 
    tline = fgetl(fileID);
    while ischar(tline)	
        if ~isempty(tline)
            if ~startsWith(tline(1), '#')
                % process lines / type and pgid should not be overwritten
                % by information from the config file
                if contains(tline,'runtype')
                    if isempty(type)
                        type = string(erase(tline,["runtype"," ","=",""""]));
                    end
                elseif contains(tline,'pgid')
                    if isempty(pgid)
                        pgid = str2num(erase(tline,["pgid"," ","="]));
                    end
                elseif contains(tline,'walltime=')
                    walltime = seconds(duration(erase(tline,["walltime"," ","="])));
                elseif contains(tline,'walltime_remaining=')
                    walltime_remaining = seconds(duration(erase(tline,["walltime_remaining"," ","="])));
                elseif contains(tline,'runtable=')
                    runtable = UserVar.home+"/"+string(erase(tline,["runtable"," ","=",""""]));
                elseif contains(tline,'idrange=')
                    idrange = str2double(split(string(erase(tline,["idrange"," ","=","""","[","]"])),":"));
                end
            end
        end
        % Now read the next line.
        tline = fgetl(fileID);
    end
    % All done reading all lines, so close the file.
    fclose(fileID);
else
    %% No config file specified, define default values
    walltime = 31557600; % set to some large number
    walltime_remaining = walltime;
    runtable = UserVar.home+"/RunTable_"+UserVar.hostname+".csv";
    idrange = [1000 1999];
end

%% check that all inputs are now available
if isempty(type) || isempty(walltime) || isempty(pgid) || isempty(walltime_remaining) || isempty(runtable)
    error("One of the following mandatory input variables is empty: type ("+string(type)+"), "+ ...
    "walltime ("+string(walltime)+"), walltime_remaining ("+string(walltime_remaining)+"), pgid ("+string(pgid)+"), ",...
    "runtable ("+string(runtable)+"), idrange ("+string(idrange)+").");
else
    UserVar.type = type;
    UserVar.pgid = pgid;
    UserVar.walltime = walltime; 
    UserVar.walltime_remaining = walltime_remaining;
    UserVar.Table = runtable;
    UserVar.idrange = idrange;
end

%% read run table
RunTable = ANT_ReadWritetable(UserVar,[],'read');

if ~isempty(RunTable)
    Iexisting = find(RunTable{:,'ExpID'}~=0 & RunTable{:,'Error'}==0);
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
        % indsnr = RunTable{ind,'Submitted'}~=0 & RunTable{ind,'Running'}==0;

        % check if not submitted, not running, not finished
        indnsnr = RunTable{ind,'Submitted'}==0 & RunTable{ind,'Running'}==0 & RunTable{ind,'Finished'}==0;
    
        % something wrong?
        % if indsnr
        %     fprintf(fid,"   ...ANT_UaWrapper: ExpID %s has been submitted, but corresponding jobID has not " + ...
        %         "been found. Either something went wrong, or the run has " + ...
        %         "finished. Check log files for errors.\n",string(RunTable{indsnr,'ExpID'}));
        %     error('');
        % end

        
        if indnsnr

            fprintf(fid,'============================\n');
            fprintf(fid,string(datetime("now"))+"\n");    
            fprintf(fid,"> ANT_UaWrapper: ExpID %s has been not yet been submitted. Let's check if a " + ...
                "restart is required...",string(RunTable{ind,'ExpID'}));

            % initialze some variables
            UserVar.Finished = 0;
            UserVar.Breakout = 0;

            % new run or restart?
            if RunTable{ind,'Restart'}==1
                UserVar.Restart = 1;
                fprintf(fid,"yes.\n");
                fprintf(fid,"> ANT_UaWrapper: Restarting ExpID %s...\n",string(RunTable{ind,'ExpID'}));
            else
                UserVar.Restart = 0;
                fprintf(fid,"no.\n");
            end

            % now gather run info and launch job
            if type=="Diagnostic"

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
                            fprintf(UserVar.fid,"> ANT_UaWrapper: Doing %s iterations or as many as walltime allows.\n",string(UserVar.TargetIterations));
    
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
                        if UserVar.Restart
                            fprintf(UserVar.fid,"> ANT_UaWrapper: %s: Breaking out of inverse cycle %s due to walltime constraints. Done %s iterations out of %s.\n",...
                            UserVar.Experiment,string(UserVar.Inverse.Cycle),string(UserVar.Inverse.IterationsDone),string(it_tmp(end)));
                        else
                            fprintf(UserVar.fid,"> ANT_UaWrapper: %s: Breaking out of inverse cycle %s. Done %s iterations out of %s.\n",...
                            UserVar.Experiment,string(UserVar.Inverse.Cycle),string(UserVar.Inverse.IterationsDone),string(it_tmp(end)));
                        end
    
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
                        fprintf(UserVar.fid,"> ANT_UaWrapper: %s: End spinup cycle %s.\n",UserVar.Experiment,string(UserVar.Spinup.Cycle));
    
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

if something_submitted % something has been submitted and run has finished -> quit matlab
    quit;
else
    RunTable = ANT_ReadWritetable(UserVar,[],'read');
    if ~isempty(RunTable)
        Inew = find(RunTable{:,'ExpID'}==0);
    else
        fprintf(fid,"Empty RunTable - stop job and quit./n");
        quit;
    end
end

% deal with new jobs if all existing jobs have been dealt with
if ~isempty(Inew)

    ind = Inew(1);

    % generate unique ExpID and save to run table
    existingID = RunTable{:,"ExpID"};
    ExpID = 0;
    while ismember(ExpID,existingID) 
        ExpID = randi([UserVar.idrange(1) UserVar.idrange(2)]);
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

    fprintf(fid,'============================\n');
    fprintf(fid,string(datetime("now"))+"\n");    
    
    % now gather run info
    if type=="Diagnostic"

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
                    fprintf(UserVar.fid,"> ANT_UaWrapper: Doing %s iterations or as many as walltime allows.\n",string(UserVar.TargetIterations));
    
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
                if UserVar.Restart
                    fprintf(UserVar.fid,"> ANT_UaWrapper: %s: Breaking out of inverse cycle %s due to walltime constraints. Done %s iterations out of %s.\n",...
                    UserVar.Experiment,string(UserVar.Inverse.Cycle),string(UserVar.Inverse.IterationsDone),string(it_tmp(end)));
                else
                    fprintf(UserVar.fid,"> ANT_UaWrapper: %s: Breaking out of inverse cycle %s. Done %s iterations out of %s.\n",...
                    UserVar.Experiment,string(UserVar.Inverse.Cycle),string(UserVar.Inverse.IterationsDone),string(it_tmp(end)));
                end
    
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
                fprintf(UserVar.fid,"> ANT_UaWrapper: %s: Breaking out of spinup cycle %s.\n",UserVar.Experiment,string(UserVar.Spinup.Cycle));
    
            end

        end
    end

    ANT_CleanUp(UserVar);

    diary off

else

    fprintf(fid,"   ...ANT_UaWrapper: Nothing to do. Try again later.\n");
    quit;

end

