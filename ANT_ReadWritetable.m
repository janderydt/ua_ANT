function RunTable = ANT_ReadWritetable(UserVar,TableToReadOrWrite,NewTable,mode)

%% Inputs:
% UserVar: needed for UserVar.home, UserVar.ExpID, UserVar.type
% TableToReadOrWrite: full <path/name> of table to read/write
% NewTable: table to write to TableToReadOrWrite - only required when
% mode=='write'
% mode: 'read' or 'write'

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

        %% Try to avoid more than 1 process is writing to the table at the same
        %% time. This is not a very robust routine, and it is generally bad
        %% practice to allow multiple processes to write to the same file.
        %% 
        isactive = 1;

        while isactive
        
            pause(1);
            isactive = isfile(strip(TableToReadOrWrite,".csv")+"_rw_active");
        
        end
        
        fid = fopen( strip(TableToReadOrWrite,".csv")+"_rw_active", 'wt' );
        fclose(fid);
        %%
      
        if isfield(UserVar,'ExpID')
            if ~isempty(UserVar.ExpID)
                % check if table already exists
                if ~isfile(TableToReadOrWrite) % if not, simply write

                    writetable(NewTable,TableToReadOrWrite);

                else % if table already exists, only replace relevant row

                    % read existing version of the table
                    opts = setoptions(UserVar);
	                RunTable_old = readtable(TableToReadOrWrite,opts);
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

            else

                writetable(NewTable,TableToReadOrWrite);

            end

        end


        %%
        delete(strip(TableToReadOrWrite,".csv")+"_rw_active");
        %%

end


end


function opts = setoptions(UserVar,TableToReadOrWrite)

opts = detectImportOptions(TableToReadOrWrite);

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
