function UserVar=DefineOutputs(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo)

%%
v2struct(F);
time=CtrlVar.time;
plots='';

if UserVar.SpinupCycle
    
    % save data in files with running names
    % check if folder 'ResultsFiles' exists, if not create   
    if exist(fullfile(cd,UserVar.UaOutputDirectory),'dir')~=7
        mkdir(UserVar.UaOutputDirectory) ;
    end
            
    FileName=sprintf('%s/%s_SpinupCycle%s_%07i.mat',...
        UserVar.UaOutputDirectory,CtrlVar.Experiment,string(UserVar.Spinup.Cycle),round(time));

    % save at y=0:1:9,10:10:100,1000:1000:end
    if round(time)==0
        fprintf(' Saving data in %s \n',FileName)
        MUA.Deriv=[]; MUA.DetJ=[]; MUA.M=[]; MUA.dm=[]; % save some space
        save(FileName,'UserVar','CtrlVar','MUA','F');
    elseif any(ismember([1:9,10:10:100,200:100:1000,2000:1000:100000],round(time)))
        if CtrlVar.AdaptMesh == 0
            fprintf(' Saving data in %s \n',FileName)
            save(FileName,'UserVar','CtrlVar','F'); % only need to save MUA once because no remeshing
        else
            fprintf(' Saving data in %s \n',FileName)
            MUA.Deriv=[]; MUA.DetJ=[]; MUA.M=[]; MUA.dm=[]; % save some space
            save(FileName,'UserVar','CtrlVar','MUA','F');
        end
    end
            
end

if contains(plots,'-plot-')
    
    figsWidth=1000 ; figHeights=300;
    GLgeo=[]; xGL=[] ; yGL=[];
    %%
    
    FindOrCreateFigure("FourPlots") ; % ,[50 50 figsWidth 3*figHeights]) ;

    subplot(4,1,1)
    PlotMeshScalarVariable(CtrlVar,MUA,F.s); title(sprintf('s at t=%g',time))
    hold on    
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL);
    %Plot_sbB(CtrlVar,MUA,s,b,B) ; title(sprintf('time=%g',time))
    
    
    subplot(4,1,2)
    QuiverColorGHG(MUA.coordinates(:,1),MUA.coordinates(:,2),F.ub,F.vb,CtrlVar);
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL);
    hold off
    
    subplot(4,1,3)
    PlotMeshScalarVariable(CtrlVar,MUA,F.dhdt);   title(sprintf('dhdt at t=%g',time))
    hold on
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL);
    
    subplot(4,1,4)
    PlotMeshScalarVariable(CtrlVar,MUA,ab);   title(sprintf('ab at t=%g',time))
    hold on
    
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL);
    hold off
    
    
    x=MUA.coordinates(:,1);
    y=MUA.coordinates(:,2);
    
    Fb=scatteredInterpolant(x,y,b);
    Fs=Fb ; Fs.Values=s;
    
    xProfile=min(x):1000:max(x);
    
    yCentre=40e3+xProfile*0;
    sProfile=Fs(xProfile,yCentre);
    bProfile=Fb(xProfile,yCentre);
    
    BProfile=MismBed(xProfile,yCentre);
    
        
    FindOrCreateFigure("Profile") ; 
    plot(xProfile/1000,sProfile,'b')
    hold on
    plot(xProfile/1000,bProfile,'b')
    plot(xProfile/1000,BProfile,'k')
    title(sprintf('t=%g',time))
    hold off
    
    
    FindOrCreateFigure("Mesh and grounding line") ; 
    PlotMuaMesh(CtrlVar,MUA);
    hold on 
    
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r','LineWidth',2);
    title(sprintf('t=%g',time))
    hold off
    
    drawnow
    %%
end

%% Deal with inverse runs
if strcmp(CtrlVar.DefineOutputsInfostring,'Start of inverse run')
    UserVar.Inverse.Iter0 =  RunInfo.Inverse.Iterations(end);
    fprintf(CtrlVar.fidlog,'(Re)Starting cycle %s of inverse run with %s iterations. Number of iterations already done: %s.\n',num2str(UserVar.Inverse.Cycle),num2str(UserVar.TargetIterations),num2str(UserVar.Inverse.Iter0));
end

if strcmp(CtrlVar.DefineOutputsInfostring,'End of Inverse Run')
    % calculate total number of iterations done over all inverse cycles
    IterationsDoneInThisRun = RunInfo.Inverse.Iterations(end)-UserVar.Inverse.Iter0;
    UserVar.Inverse.IterationsDone = UserVar.Inverse.IterationsDone + IterationsDoneInThisRun;
	if IterationsDoneInThisRun ~= CtrlVar.Inverse.Iterations
        if RunInfo.Inverse.stoppedduetowalltime == 1
            fprintf(CtrlVar.fidlog,['Simulation stopped due to walltime constraints. Done %s iterations instead of %s. ',...
        		'Writing restart file.\n'],num2str(IterationsDoneInThisRun),num2str(CtrlVar.Inverse.Iterations));         
            UserVar.Restart = 1;
            UserVar.Error = 0;            
        else
    		fprintf(CtrlVar.fidlog,['Simulation did not reach expected number of iterations. Done %s instead of %s. ',...
        		'Writing restart file and breaking out.\n'],num2str(IterationsDoneInThisRun),num2str(CtrlVar.Inverse.Iterations));
            RunTable_exp=ANT_ReadWritetable(UserVar,UserVar.runtable_exp,[],'read');
            ind = find(RunTable_exp{:,'ExpID'}==UserVar.ExpID);
            RunTable_exp{ind,'Comments'}={['Simulation did not reach expected number of iterations. Done ',...
                num2str(IterationsDoneInThisRun),' instead of ',num2str(CtrlVar.Inverse.Iterations),'.']};
            [~]=ANT_ReadWritetable(UserVar,UserVar.runtable_exp,RunTable_exp,'write');
            UserVar.Restart = 0;
            UserVar.Error = 1;
        end     
        
    else
        fprintf(CtrlVar.fidlog,'Simulation reached expected number of %s iterations.\n',num2str(IterationsDoneInThisRun));
        UserVar.Restart = 0;
        UserVar.Error = 0;
    end
    UserVar.Finished = 0;
    UserVar.Breakout = 1;

    WriteAdjointRestartFile(UserVar,CtrlVar,MUA,BCs,F,F.GF,l,RunInfo,InvStartValues,Priors,Meas,BCsAdjoint,InvFinalValues);
    NameOfRestartOutputFile = erase(CtrlVar.Inverse.NameOfRestartOutputFile,".mat")+"_InverseCycle"+string(UserVar.Inverse.Cycle)+".mat";
    copyfile(CtrlVar.Inverse.NameOfRestartOutputFile,NameOfRestartOutputFile);
    NameOfRestartOutputFile = erase(CtrlVar.Inverse.NameOfRestartOutputFile,".mat")+"_It"+string(UserVar.Inverse.IterationsDone)+".mat";
    copyfile(CtrlVar.Inverse.NameOfRestartOutputFile,NameOfRestartOutputFile);
end

%% Deal with spinup simulations
if strcmp(CtrlVar.DefineOutputsInfostring,'Last call') % the string "last call" is only set for a transient simulation, not for an inverse simulation
    % calculate total number of years over all spinup cycles 
    UserVar.Spinup.YearsDone = CtrlVar.time;
    UserVar.Error = 0;
    UserVar.Finished = 0;
    UserVar.Breakout = 1;
    if UserVar.Spinup.stoppedduetowalltime == 1
        fprintf(CtrlVar.fidlog,['Simulation stopped due to walltime constraints. Done %s years instead of %s. ',...
        		'Writing restart file.\n'],num2str(UserVar.Spinup.YearsDone),num2str(CtrlVar.TotalTime));         
        UserVar.Restart = 1;         
    else
        if isapprox(UserVar.Spinup.YearsDone,CtrlVar.TotalTime)
            fprintf(CtrlVar.fidlog,'Simulation reached expected number of %s years.\n',num2str(CtrlVar.TotalTime));            
            UserVar.Restart = 0;
        else
            % simulation finished before the walltime without error, but 
            % did not reach expected number of years. One reason can be
            % that the timestep became too small
            RunTable_exp=ANT_ReadWritetable(UserVar,UserVar.runtable_exp,[],'read');
            ind = find(RunTable_exp{:,'ExpID'}==UserVar.ExpID);
            if CtrlVar.dt<=CtrlVar.dtmin   
                RunTable_exp{ind,'Comments'}={['Simulation did not reach expected number of ',...
                    num2str(CtrlVar.TotalTime),'. Timestep too small.']};          
            else
                RunTable_exp{ind,'Comments'}={['Simulation did not reach expected number of ',...
                    num2str(CtrlVar.TotalTime),'. Check log file for reason.']};
            end  
            [~]=ANT_ReadWritetable(UserVar,UserVar.runtable_exp,RunTable_exp,'write');
            UserVar.Restart = 0;
            UserVar.Error = 1;
        end
    end
    WriteForwardRunRestartFile(UserVar,CtrlVar,MUA,BCs,F,GF,l,RunInfo);
    NameOfRestartOutputFile = erase(CtrlVar.NameOfRestartFiletoWrite,".mat")+"_SpinupCycle"+string(UserVar.Spinup.Cycle)+".mat";
    copyfile(CtrlVar.NameOfRestartFiletoWrite,NameOfRestartOutputFile);
    NameOfRestartOutputFile = erase(CtrlVar.NameOfRestartFiletoWrite,".mat")+"_Yrs"+strrep(string(CtrlVar.time),".","k")+".mat";
    copyfile(CtrlVar.NameOfRestartFiletoWrite,NameOfRestartOutputFile);
end

%% Deal with sequence of inverse and diagnostic runs
if contains(CtrlVar.DefineOutputsInfostring,{'Last call','End of Inverse Run'})
    years_tmp = cumsum(UserVar.Spinup.Years);
    it_tmp = cumsum(UserVar.Inverse.Iterations);
    if UserVar.Inverse.IterationsDone == it_tmp(end) && UserVar.Spinup.YearsDone == years_tmp(end)
        UserVar.Finished = 1;
        UserVar.Restart = 0;
        UserVar.Error = 0;
        fprintf(CtrlVar.fidlog,'Simulation completed %s inverse iterations in %s cycles and %s spinup years in %s cycles. Breaking out.\n',...
            num2str(UserVar.Inverse.IterationsDone),num2str(UserVar.Inverse.Cycle),num2str(UserVar.Spinup.YearsDone),num2str(UserVar.Spinup.Cycle));
    end      
end

end
