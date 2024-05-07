function RunTable = ANT_ReadWritetable(UserVar,RunTable,mode)

switch mode 
    case 'read'
        if exist(UserVar.Table,'file')
            opts = detectImportOptions(UserVar.Table);
            % set some variable types    
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
            
            %I=find(contains(opts.VariableNames,'SubmissionTime'));
            %opts.VariableTypes(I)={'char'};
            %opts = setvaropts(opts,'SubmissionTime','InputFormat','dd/MM/yyyy HH:mm:ss'); 

            %I=find(contains(opts.VariableNames,'ErrorTime'));
            %opts.VariableTypes(I)={'char'};
            %opts = setvaropts(opts,'ErrorTime','InputFormat','dd/MM/yyyy HH:mm:ss'); 

            %I=find(contains(opts.VariableNames,'FinishedTime'));
            %opts.VariableTypes(I)={'char'};
            %opts = setvaropts(opts,'FinishedTime','InputFormat','dd/MM/yyyy HH:mm:ss'); 

            % read table
            RunTable=readtable(UserVar.Table,opts); 
        else    
            error("Runtable does not exist"); 
        end
    case 'write'
        writetable(RunTable,UserVar.Table);
end
