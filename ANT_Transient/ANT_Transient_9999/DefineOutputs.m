function UserVar=DefineOutputs(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo)
%%
% This routine is called during the run and can be used for saving and/or plotting data.
% 
%   UserVar=DefineOutputs(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo)
%
% Write your own version of this routine and put it in you local run directory.
% 
%
%   This is the m-file you use to define/plot your results.
%
%   You will find all the outputs in the variable F
%
%   The variable F is a structure, and has various fields.
%
%   For example:
%
%   F.s             The upper glacier surface%
%   F.b             The lower glacier surface
%   F.B             The bedrock
%   F.rho           The ice density
%   F.C             Basal slipperiness, i.e. the variable C in the basal sliding law
%   F.AGlen         The rate factor, i.e. the variable A in Glen's flow law
%
%   F.ub            basal velocity in x-direction
%   F.vb            basal velocity in y-direction
%
%   All these variables are nodal variables, i.e. these are the corresponding values at the nodes of the computational domain.
%
%   You find informaton about the computational domain in the variable MUA
%
%   For example, the x and y coordinates of the nodes are in the nx2 array MUA.coordinates, where n is the number of nodes.
%
%   MUA.coordinates(:,1)    are the nodal x coordinates
%   MUA.coordinates(:,y)    are the nodal y coordinates
%
%
%   BCs             Structure with all boundary conditions
%   l               Lagrange parameters related to the enforcement of boundary
%                   conditions.
%   GF              Grounding floating mask for nodes and elements.
%
%   Note: If preferred to work directly with the variables rather than the respective fields of the structure F, thenF can easily be
%   converted into variables using v2struc.
%
%
%
%   Note:  For each call to this m-File, the variable
%
%       CtrlVar.DefineOutputsInfostring
%
%   gives you information about different stages of the run (start, middle
%   part, end, etc.).
%
%   So for example, when Ua calls this m-File for the last time during the
%   run, the variable has the value
%
%     CtrlVar.DefineOutputsInfostring="Last call"
%
%
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
    
    if UserVar.stoppedduetowalltime == 1
        fprintf(CtrlVar.fidlog,['Simulation stopped due to walltime constraints. Done %s years instead of %s. ',...
        		'Writing restart file.\n'],num2str(UserVar.YearsDone),num2str(CtrlVar.TotalTime));         
        UserVar.Restart = 1;   
        UserVar.Finished = 0;     
    else
        fprintf(CtrlVar.fidlog,'Simulation reached expected number of %s years.\n',num2str(CtrlVar.TotalTime));
        UserVar.Restart = 0; 
        UserVar.Finished = 1;
    end
    UserVar.Error = 0;    
    UserVar.Breakout = 1;
    WriteForwardRunRestartFile(UserVar,CtrlVar,MUA,BCs,F,GF,l,RunInfo);  
    NameOfRestartOutputFile = erase(CtrlVar.NameOfRestartFiletoWrite,".mat")+"_Yrs"+strrep(string(CtrlVar.time),".","k")+".mat";
    copyfile(CtrlVar.NameOfRestartFiletoWrite,NameOfRestartOutputFile);
end


end
