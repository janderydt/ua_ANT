function ManuallyAdjustRunTable

% It might be necessary to adjust an existing runtable during run-time, 
% based on user-defined criteria. Rather than doing this in Excell, it is
% much easier to write a little script. One example I have used this for,
% is to let inverse simulations that have stagnated continue to the spin-up
% phase, even though they have not reached the required number of inverse
% steps.

addpath("/mnt/md0/Ua/cases/ANT/");
addpath(getenv("froot_tools"));

UserVar.home = "/mnt/md0/Ua/cases/ANT/";
UserVar.type = "Inverse";
UserVar.domain = "AMUND";

RunTables_to_adjust = UserVar.home+"ANT_"+UserVar.type+"/RunTable_ARCHER2_17-02-2025_"+string([25])+".csv";

for tt=1:numel(RunTables_to_adjust)

    %% read run table
    RunTable = ANT_ReadWritetable(UserVar,RunTables_to_adjust(tt),[],'read');

    % error flag:
    errorflag = RunTable{:,'Error'};
    % target iterations in first inversion:
    targetiterations_full = string(cell2mat(RunTable{:,'InverseIterations'}));
    targetiterations = [];
    for ii=1:numel(targetiterations_full)
        tmp = strsplit(targetiterations_full(ii),"+");
        targetiterations(ii) = double(tmp(1));
    end
    % inversion steps completed so far:
    iterationsdone = RunTable{:,'InverseIterationsDone'};
    
    % find runs that flagged error, and for which less than the target
    % iterations were completed
    Ind = find(errorflag==1 & iterationsdone(:)<targetiterations(:));

    % for those simulations, set errorflag to False, and replace target
    % iterations to number of completed interations

    targetiterations_full_new = targetiterations_full;
    errorflag_new = errorflag;
    for ii=1:numel(Ind)
        targetiterations_full_new(Ind(ii)) = strrep(targetiterations_full_new(Ind(ii)),...
            string(targetiterations(Ind(ii))),string(iterationsdone(Ind(ii))));
        errorflag_new(Ind(ii)) = 0;
        RunTable{Ind(ii),'InverseIterations'}{:} = char(targetiterations_full_new(Ind(ii)));
        RunTable{Ind(ii),'Error'} = errorflag_new(Ind(ii));      
        RunTable{Ind(ii),'Submitted'} = 0;      
    end

    % save adjusted RunTable
    ANT_ReadWritetable(UserVar,RunTables_to_adjust(tt),RunTable,'write');
end