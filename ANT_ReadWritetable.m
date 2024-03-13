function RunTable = ANT_ReadWritetable(UserVar,RunTable,mode)

switch mode 
    case 'read'
        if exist(UserVar.Table,'file')
            opts = detectImportOptions(UserVar.Table);
            % set some variable types    
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
            % read table
            RunTable=readtable(UserVar.Table,opts);  
        else    
            error("Runtable does not exist"); 
        end
    case 'write'
        writetable(RunTable,UserVar.Table);
end
