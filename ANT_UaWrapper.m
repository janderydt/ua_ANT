function ANT_UaWrapper(ua_config,pgid,type,row_number,expid_new)

%% INPUTS
% > ua_config: link to *txt file with config details such as the run table,
%   run type (Inverse, Diagnostic, Transient), total walltime and remaining 
%   walltime / format: string
% > pgid: on linux workstations, each job has a unique pgid, which can be
%   added as an argument here / format: integer
% > type: when the ua_config file is not specified, provide the run type
%   (Inverse, Diagnostic, Transient) as an argument / format: string
% > row_number is the row number in the RunTable of the experiment to 
%   submit / format: integer
% > expid: when row_number is specified, you also need to provide a unique 
%  experiment id as an argument / format: integer
arguments
    ua_config (1,1) string = ''
    pgid (1,1) string = ''
    type (1,1) string = ''
    row_number (1,1) string = ''
    expid_new (1,1) string = ''
end

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

UserVar.home = pwd+"/";

%% deal with inputs
if nargin == 0
    error("You need to either provide a config file or pgid and type. Instead got none.");
elseif nargin == 4
    error("You have specified a row number, but also need to provide a unique experiment id for that row.")
end
if nargin > 1
    pgid = str2double(pgid);
end
if nargin > 3
    row_number = round(double(row_number));
    expid_new = round(double(expid_new));
end
    
if ~isempty(ua_config)
    % read inputs from config file
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
                    if strlength(type)==0
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
    % No config file specified, define default values
    walltime = 31557600; % set to some large number
    walltime_remaining = walltime;
    runtable = UserVar.home+"/RunTable_"+UserVar.hostname+".csv";
    idrange = [1 999];
end

% check that all inputs are now available
if isempty(type) || isempty(walltime) || isempty(pgid) || isempty(walltime_remaining) || isempty(runtable)
    error("One of the following mandatory input variables is empty: type ("+string(type)+"), "+ ...
    "walltime ("+string(walltime)+"), walltime_remaining ("+string(walltime_remaining)+"), pgid ("+string(pgid)+"), ",...
    "runtable ("+string(runtable)+"), idrange ("+string(idrange)+").");
else
    UserVar.type = type;
    UserVar.pgid = pgid;
    UserVar.walltime = walltime; 
    UserVar.walltime_remaining = walltime_remaining;
    UserVar.runtable_global = runtable;
    UserVar.idrange = idrange;
end

%% initialize global log file
logfile = UserVar.home+"/jobs_master_"+UserVar.hostname+".log";
fid = fopen(logfile,'a+');
UserVar.fid_masterlog = fid;

%% read run table
RunTable = ANT_ReadWritetable(UserVar,UserVar.runtable_global,[],'read');

if ~isempty(RunTable)
    Iexisting = find(RunTable{:,'ExpID'}~=0 & RunTable{:,'Error'}==0);
    if ~isempty(Iexisting) & ~isempty(row_number)
        % check if row_number is an existing simulation
        if ismember(row_number,Iexisting)
            Iexisting = row_number;
        else
            Iexisting = [];
        end
    end
end

%% launch jobs
% initialize some variables
UserVar.configfile = ua_config;
UserVar.casefolder = UserVar.home + "/cases/";
UserVar.datafolder = UserVar.home + "/../ANT_Data/";
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

        % initialize experiment log file
        logfile = UserVar.casefolder+"/"+string(UserVar.Experiment)+"/"+string(UserVar.Experiment)+".log";
        fid = fopen(logfile,'a+');
        UserVar.fid_experimentlog = fid;

        % check if not submitted, not running, not finished
        indnsnr = RunTable{ind,'Submitted'}==0 & RunTable{ind,'Running'}==0 & RunTable{ind,'Finished'}==0;
        
        if indnsnr

            fprintf(UserVar.fid_experimentlog,'============================\n');
            fprintf(UserVar.fid_experimentlog,string(datetime("now"))+"\n");    
            fprintf(UserVar.fid_experimentlog,"> ANT_UaWrapper: ExpID %s has been not yet been submitted. Let's check if a " + ...
                "restart is required...",string(RunTable{ind,'ExpID'}));

            % initialze some variables
            UserVar.Finished = 0;
            UserVar.Breakout = 0;
            UserVar.runtable_exp = UserVar.casefolder + "/" + UserVar.Experiment + "/RunTable_" + UserVar.Experiment + ".csv";

            % new run or restart?
            if RunTable{ind,'Restart'}==1
                UserVar.Restart = 1;
                fprintf(UserVar.fid_experimentlog,"yes.\n");
                fprintf(UserVar.fid_experimentlog,"> ANT_UaWrapper: Restarting ExpID %s...\n",string(RunTable{ind,'ExpID'}));
            else
                UserVar.Restart = 0;
                fprintf(UserVar.fid_experimentlog,"no.\n");
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
                    RunTable=ANT_ReadWritetable(UserVar,UserVar.runtable_global,[],'read');
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
                                nit = 100000;
                            else
                                nit = 5000;
                            end
        
                            UserVar.TargetIterations = min(nit,it_tmp(UserVar.Inverse.Cycle)-UserVar.Inverse.IterationsDone);
                            fprintf(UserVar.fid_experimentlog,"> ANT_UaWrapper: Doing %s iterations or as many as walltime allows.\n",string(UserVar.TargetIterations));
    
                            UserVar = ANT_UaJob(RunTable,ind,UserVar,pgid);
    
                            %adjust Individual Runtable
                            RunTable=ANT_ReadWritetable(UserVar,UserVar.runtable_exp,[],'read');
                            ind = find(RunTable{:,'ExpID'}(:) == UserVar.ExpID);
                            RunTable{ind,"InverseIterationsDone"} = UserVar.Inverse.IterationsDone;   
                            RunTable{ind,"Restart"} = UserVar.Restart;
                            [~] = ANT_ReadWritetable(UserVar,UserVar.runtable_exp,RunTable,'write');
                
                        end   

                        it_tmp = cumsum(UserVar.Inverse.Iterations);

                        fprintf(UserVar.fid_experimentlog,'============================\n');
                        fprintf(UserVar.fid_experimentlog,string(datetime("now"))+"\n");                     
                        if UserVar.Restart
                            fprintf(UserVar.fid_experimentlog,"> ANT_UaWrapper: %s: Breaking out of inverse cycle %s due to walltime constraints. Done %s iterations out of %s.\n",...
                            UserVar.Experiment,string(UserVar.Inverse.Cycle),string(UserVar.Inverse.IterationsDone),string(it_tmp(end)));
                        else
                            fprintf(UserVar.fid_experimentlog,"> ANT_UaWrapper: %s: Breaking out of inverse cycle %s. Done %s iterations out of %s.\n",...
                            UserVar.Experiment,string(UserVar.Inverse.Cycle),string(UserVar.Inverse.IterationsDone),string(it_tmp(end)));
                        end
    
                    %% Spinup cycle
                    elseif UserVar.SpinupCycle
    
                        UserVar = ANT_UaJob(RunTable,ind,UserVar,pgid);
    
                        %adjust Runtable
                        RunTable=ANT_ReadWritetable(UserVar,UserVar.runtable_exp,[],'read');
                        ind = find(RunTable{:,'ExpID'}(:) == UserVar.ExpID);
                        RunTable{ind,"SpinupYearsDone"} = UserVar.Spinup.YearsDone;   

                        [~] = ANT_ReadWritetable(UserVar,UserVar.runtable_exp,RunTable,'write');

                        fprintf(UserVar.fid_experimentlog,'============================\n');
                        fprintf(UserVar.fid_experimentlog,string(datetime("now"))+"\n");                        
                        fprintf(UserVar.fid_experimentlog,"> ANT_UaWrapper: %s: End spinup cycle %s.\n",UserVar.Experiment,string(UserVar.Spinup.Cycle));
    
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
    RunTable = ANT_ReadWritetable(UserVar,UserVar.runtable_global,[],'read');
    if ~isempty(RunTable)
        Inew = find(RunTable{:,'ExpID'}==0);
        if ~isempty(Inew) & ~isempty(row_number)
            % check if row_number is a new simulation
            if ismember(row_number,Inew)
                Inew = row_number;
                ExpID = expid_new;
            else
                Inew = [];
                ExpID = 0;
            end
        end
    else
        fprintf(UserVar.fid_masterlog,"Empty RunTable - stop job and quit./n");
        quit;
    end
end

% deal with new jobs if all existing jobs have been dealt with
if ~isempty(Inew)

    ind = Inew(1);

    % if expid hasn't been specified as input then generate unique ExpID
    % and write to global runtable
    if ExpID == 0
        existingID = RunTable{:,"ExpID"};
        while ismember(ExpID,existingID) 
            ExpID = randi([UserVar.idrange(1) UserVar.idrange(2)]);
        end   
        RunTable{ind,"ExpID"} = ExpID;
        RunTable{ind,'pgid'} = pgid;
        [~]=ANT_ReadWritetable(UserVar,UserVar.runtable_global,RunTable,'write');
    else
        RunTable{ind,"ExpID"} = ExpID;
        RunTable{ind,'pgid'} = pgid;
    end

    % initialize some UserVars
    UserVar.ExpID = ExpID;
    UserVar.Finished = 0;
    UserVar.Restart = 0;
    UserVar.Breakout = 0;
    UserVar.Domain = RunTable{ind,'Domain'};
    UserVar.Experiment = [char(UserVar.Domain),'_',char(type),'_',num2str(ExpID)];
    UserVar.runtable_exp = UserVar.casefolder + "/" + UserVar.Experiment + "/RunTable_" + UserVar.Experiment + ".csv";

    % make copy of master folder for new experiment
    % if new folder already exists: rename first
    sourcefolder = [char(UserVar.home),'/ANT_',char(type),'_9999/'];
    newfolder = [char(UserVar.casefolder),'/',char(UserVar.Domain),'_',char(type),'_',num2str(ExpID)];
    if exist(newfolder,"dir") == 7
        movefile(newfolder,[newfolder,'_old/']);
    else
        copyfile(sourcefolder,[newfolder,'/']); 
    end
    % rename RunTable file
    movefile(newfolder+"/RunTable_ANT_"+UserVar.type+"_9999.csv",UserVar.runtable_exp);

    % initialize experiment log file
    logfile = UserVar.casefolder+"/"+string(UserVar.Experiment)+"/"+string(UserVar.Experiment)+".log";
    fid = fopen(logfile,'a+');
    UserVar.fid_experimentlog = fid;

    fprintf(UserVar.fid_experimentlog,'============================\n');
    fprintf(UserVar.fid_experimentlog,string(datetime("now"))+"\n");    
    
    % now gather run info
    if type=="Diagnostic"

        % initialize User variables
        UserVar = ANT_GetUserVar_Diagnostic(RunTable,ind,UserVar);

        UserVar.Restart = 0;
    
        UserVar = ANT_UaJob(RunTable,ind,UserVar,pgid);

    elseif type=="Inverse"

        while ~UserVar.Breakout
            
            % read Runtable again in case any changes were made by other
            % processes
            %RunTable=ANT_ReadWritetable(UserVar,UserVar.runtable_global,[],'read');
            %ind = find(RunTable{:,'ExpID'}(:) == UserVar.ExpID);
               
            % initialize User variables
            UserVar = ANT_GetUserVar_Inverse(RunTable,ind,UserVar);
    
            % cumulative sum of number of iterations at the end of each
            % inverse cycle
            it_tmp = cumsum(UserVar.Inverse.Iterations);
    
            %% Inverse cycle
            if UserVar.InverseCycle
            
                while (UserVar.Inverse.IterationsDone < it_tmp(UserVar.Inverse.Cycle) && ~UserVar.Breakout)
    
                    if contains(UserVar.hostname,"ARCHER2")
                        nit = 100000;
                    else
                        nit = 5000;
                    end

                    UserVar.TargetIterations = min(nit,it_tmp(UserVar.Inverse.Cycle)-UserVar.Inverse.IterationsDone);
                    fprintf(UserVar.fid_experimentlog,"> ANT_UaWrapper: Doing %s iterations or as many as walltime allows.\n",string(UserVar.TargetIterations));
    
                    UserVar = ANT_UaJob(RunTable,ind,UserVar,pgid);
    
                    %adjust Runtable
                    RunTable=ANT_ReadWritetable(UserVar,UserVar.runtable_exp,[],'read');
                    ind = find(RunTable{:,'ExpID'}(:) == UserVar.ExpID);
                    RunTable{ind,"InverseIterationsDone"} = UserVar.Inverse.IterationsDone;  
                    RunTable{ind,"Restart"} = UserVar.Restart;
                    [~] = ANT_ReadWritetable(UserVar,UserVar.runtable_exp,RunTable,'write');

                end     

                it_tmp = cumsum(UserVar.Inverse.Iterations);

                fprintf(UserVar.fid_experimentlog,'============================\n');
                fprintf(UserVar.fid_experimentlog,string(datetime("now"))+"\n");                     
                if UserVar.Restart
                    fprintf(UserVar.fid_experimentlog,"> ANT_UaWrapper: %s: Breaking out of inverse cycle %s due to walltime constraints. Done %s iterations out of %s.\n",...
                    UserVar.Experiment,string(UserVar.Inverse.Cycle),string(UserVar.Inverse.IterationsDone),string(it_tmp(end)));
                else
                    fprintf(UserVar.fid_experimentlog,"> ANT_UaWrapper: %s: Breaking out of inverse cycle %s. Done %s iterations out of %s.\n",...
                    UserVar.Experiment,string(UserVar.Inverse.Cycle),string(UserVar.Inverse.IterationsDone),string(it_tmp(end)));
                end
    
            %% Spinup cycle
            elseif UserVar.SpinupCycle
    
                UserVar = ANT_UaJob(RunTable,ind,UserVar,pgid);

                %adjust Runtable
                RunTable=ANT_ReadWritetable(UserVar,UserVar.runtable_exp,[],'read');
                ind = find(RunTable{:,'ExpID'}(:) == UserVar.ExpID);
                RunTable{ind,"SpinupYearsDone"} = UserVar.Spinup.YearsDone;   
                [~] = ANT_ReadWritetable(UserVar,UserVar.runtable_exp,RunTable,'write');

                fprintf(UserVar.fid_experimentlog,'============================\n');
                fprintf(UserVar.fid_experimentlog,string(datetime("now"))+"\n");                
                fprintf(UserVar.fid_experimentlog,"> ANT_UaWrapper: %s: Breaking out of spinup cycle %s.\n",UserVar.Experiment,string(UserVar.Spinup.Cycle));
    
            end

        end
    end

    ANT_CleanUp(UserVar);

    diary off

else

    fid2 = fopen( UserVar.home+"/"+string(UserVar.pgid)+"_job_submitted", 'wt' );
    fclose(fid2);
    %fprintf(UserVar.fid_masterlog,"   ...ANT_UaWrapper: Nothing to do. Try again later.\n");

end

% At the end always exit matlab to avoid rogue jobs
quit;

