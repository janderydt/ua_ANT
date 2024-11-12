function UserVar=DefineOutputs(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo)
%%

v2struct(F);

UserVar.YearsCompleted = CtrlVar.time-UserVar.StartTime_DecimalYears;

plots='-save-';

%% SAVE outputs
if contains(plots,'-save-')
    
    % save data in files with running names
    % check if folder 'ResultsFiles' exists, if not create
    
    if exist(fullfile(cd,UserVar.UaOutputDirectory),'dir')~=7
        mkdir(UserVar.UaOutputDirectory) ;
    end
       
    FileName=sprintf("%s/ResultsFile-%s-0101%s-%s.mat",...
            UserVar.UaOutputDirectory,UserVar.Experiment,...
            string(UserVar.StartYear),num2str(round(time*365.25),'%06.f'));
    fprintf(' Saving data in %s \n',FileName)
    save(FileName,'UserVar','CtrlVar','MUA','F');
    
end

%% At the END of the simulation:
if strcmp(CtrlVar.DefineOutputsInfostring,'Last call') % the string "last call" is only set for a transient simulation, not for an inverse simulation
    % calculate total number of years over all spinup cycles 

    UserVar.Error = 0;    
    UserVar.Breakout = 1;
    if UserVar.stoppedduetowalltime == 1
        fprintf(CtrlVar.fidlog,['Simulation stopped due to walltime constraints. Done %s years instead of %s. ',...
        		'Writing restart file.\n'],num2str(UserVar.YearsDone),num2str(CtrlVar.TotalTime));         
        UserVar.Restart = 1;   
        UserVar.Finished = 0;     
    else
        if isapprox(UserVar.Spinup.YearsDone,CtrlVar.TotalTime)
            fprintf(CtrlVar.fidlog,'Simulation reached expected number of %s years.\n',num2str(CtrlVar.TotalTime));            
            UserVar.Restart = 0;
            UserVar.Finished = 1;     
        else
            % simulation finished before the walltime without error, but 
            % did not reach expected number of years. One reason can be
            % that the timestep became too small
            RunTable_exp=ANT_ReadWritetable(UserVar,UserVar.runtable_exp,[],'read');
            ind = find(RunTable{:,'ExpID'}==UserVar.ExpID);
            if CtrlVar.dt<=CtrlVar.dtmin   
                RunTable_exp{ind,"Comments"}="Simulation did not reach expected number of "+...
                    num2str(num2str(CtrlVar.TotalTime))+". Timestep too small.";          
            else
                RunTable_exp{ind,"Comments"}="Simulation did not reach expected number of "+...
                    num2str(num2str(CtrlVar.TotalTime))+". Check log file for reason.";
            end  
            [~]=ANT_ReadWritetable(UserVar,UserVar.runtable_exp,RunTable_exp,'write');
            UserVar.Finished = 0;
            UserVar.Restart = 0;
            UserVar.Error = 1;
        end
    end
    WriteForwardRunRestartFile(UserVar,CtrlVar,MUA,BCs,F,GF,l,RunInfo);  
    NameOfRestartOutputFile = erase(CtrlVar.NameOfRestartFiletoWrite,".mat")+"_Yrs"+strrep(string(CtrlVar.time),".","k")+".mat";
    copyfile(CtrlVar.NameOfRestartFiletoWrite,NameOfRestartOutputFile);
end


end
