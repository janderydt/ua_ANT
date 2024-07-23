function RunTable = ANT_ReadWritetable(UserVar,TableToReadOrWrite,NewTable,mode)

%% Inputs:
% UserVar: needed for UserVar.home, UserVar.ExpID, UserVar.type
% TableToReadOrWrite: full <path/name> of table to read/write
% NewTable: table to write to TableToReadOrWrite - only required when
% mode=='write'
% mode: 'read' or 'write'

RunTable = [];

switch mode 
    case 'read'
        if exist(TableToReadOrWrite,'file')
            % set some variable types    
            opts = setoptions(UserVar,TableToReadOrWrite);

            % read table
            RunTable=readtable(TableToReadOrWrite,opts);
        else
            error("Runtable "+string(TableToReadOrWrite)+" does not exist");
        end

    case 'write'
      
        if isfield(UserVar,'ExpID')

            if ~isempty(UserVar.ExpID)

                % check if table already exists
                if ~isfile(TableToReadOrWrite) % if not, simply write

                    writetable(NewTable,TableToReadOrWrite);

                else % if table already exists, only replace relevant row

                    % read existing version of the table
                    opts = setoptions(UserVar,TableToReadOrWrite);
	                RunTable_old = readtable(TableToReadOrWrite,opts);

                    if numel(RunTable_old)==0

                        writetable(NewTable,TableToReadOrWrite);

                    else
                    
                        ind_old = find(RunTable_old{:,'ExpID'}(:) == UserVar.ExpID);

                        if ~isempty(ind_old)
    
                            % only write the line that matters for this experiment, not the
                            % whole table
                            ind_new = find(NewTable{:,'ExpID'}(:) == UserVar.ExpID);  
                            % replace 1 line and write
                            RunTable_old(ind_old,:) = NewTable(ind_new,:);
    
                            writetable(RunTable_old,TableToReadOrWrite);
    
                        else
    
                            writetable(NewTable,TableToReadOrWrite);
    
                        end
   
                    end

                end

            else

                writetable(NewTable,TableToReadOrWrite);

            end

        else

            writetable(NewTable,TableToReadOrWrite);

        end

end


end


function opts = setoptions(UserVar,TableToReadOrWrite)

opts = detectImportOptions(TableToReadOrWrite);

I=find(contains(opts.VariableNames,'SubmissionTime'));
opts.VariableTypes(I)={'char'};
I=find(contains(opts.VariableNames,'ErrorTime'));
opts.VariableTypes(I)={'char'};
I=find(contains(opts.VariableNames,'FinishedTime'));
opts.VariableTypes(I)={'char'};
% opts=setvaropts(opts,'SubmissionTime','InputFormat','dd/MM/uuuu HH:mm:ss');
% opts=setvaropts(opts,'ErrorTime','InputFormat','dd/MM/uuuu HH:mm:ss');
% opts=setvaropts(opts,'FinishedTime','InputFormat','dd/MM/uuuu HH:mm:ss');

if UserVar.type == "Inverse"

    I=find(contains(opts.VariableNames,'Domain'));
    opts.VariableTypes(I)={'char'};
    I=find(contains(opts.VariableNames,'InverseIterations'));
    opts.VariableTypes(I)={'char'};
    I=find(contains(opts.VariableNames,'InverseIterationsDone'));
    opts.VariableTypes(I)={'double'};
    I=find(contains(opts.VariableNames,'SpinupYears'));
    opts.VariableTypes(I)={'char'};
    I=find(contains(opts.VariableNames,'SpinupYearsDone'));
    opts.VariableTypes(I)={'double'};
    I=find(contains(opts.VariableNames,'Comments'));
    opts.VariableTypes(I)={'char'};
                
elseif UserVar.type == "Diagnostic"

    % none to correct

elseif UserVar.type == "Transient"

    % none to correct

end

end
