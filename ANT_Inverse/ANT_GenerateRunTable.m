function ANT_GenerateRunTable(X,GradientCalc,SlidingLaw,Enrich,UseCatalogue)

%% INPUTS:
% -> X is the runtable created by ANT_GenerateRunTable or
%    ANT_ExpandExistingRunTable
% -> GradientCalc is either 'fixpoint' or 'adjoint'
% -> SlidingLaw is any valid type of sliding law in Ua, e.g. Weertman
% -> Enrich is 0 or 1, it only affects the name of the csv file that is
%    produced

addpath('../');

MaxTableLength = 270;
UserVar.home = pwd;
DT=datetime;
DT.Format="dd-MM-yyyy";
UserVar.Table = "NewTable_"+GradientCalc+"_"+SlidingLaw+"_Enriched"+string(Enrich)+"_"+string(DT)+".csv";
UserVar.type = "Inverse";

if ~exist(UserVar.Table,"file")
    copyfile("EmptyTable.csv",UserVar.Table);
end

RunTable = ANT_ReadWritetable(UserVar,UserVar.Table,[],'read');

for ind=1:size(X,1)

    startMesh = "2000_2009_2014_2018_meshmin3000_meshmax100000_refined";

    m = round(X.m(ind)*100)/100;
    n = round(X.n(ind)*100)/100; 
    gsC = round(X.gsC(ind)*10)/10;         %gsC
    gsA = round(X.gsA(ind)*10)/10;         %gsA
    gaC = round(X.gaC(ind)*10)/10;         %gaC
    gaA = round(X.gaA(ind)*10)/10;         %gaA
    dhdt_err = round(X.dhdt_err(ind)*100)/100;

    taub = 80;  
    ub = 100;%X(ind,6);
    priorC = ub/taub^m;
    muk = 0.5;
    
    tau = 80; 
    eps = 0.0026;%X(ind,8);
    priorAGlen = eps/tau^n;
   
    % If GradientCalc=Adjoint, and SlidingLaw=Weertman and Enrich=0 then 
    % use Fixpoint inversion with nearest (m,n) pair as start value. 
    % The code will rescale C and AGlen depending on m and n.
    % If GradientCalc=Adjoint, and SlidingLaw=Umbi and/or Enrich=1 then 
    % use finished Weertman inversion with nearest (m,n,gsA,gsC,gaA,gaC) 
    % as start value. The code will rescale C and AGlen depending on m and
    % n.
    if GradientCalc == "Adjoint" && Enrich == 0 && UseCatalogue == 0

        if ind == 1
            UserVar2.type = "Inverse";
            UserVar2.home = pwd; 
            UserVar2.Table = "RunTable_sauron_FixPoint.csv";
            RunTable_FixPoint = ANT_ReadWritetable(UserVar2,UserVar2.Table,[],'read');
            % scan table for ExpID, m and n
            ExpID_FixPoint = RunTable_FixPoint{:,'ExpID'};
            m_FixPoint = RunTable_FixPoint{:,'m'};
            n_FixPoint = RunTable_FixPoint{:,'n'};
        end
        [~,Ind] = min(hypot(m_FixPoint-m,n_FixPoint-n));
        startC = ExpID_FixPoint(Ind);
        startAGlen = 0;
        iterations = "15000+1000";
        spinupyears = "1";
        invertfor = "-logC-logA-";

     elseif GradientCalc == "Adjoint" && (Enrich == 1 || UseCatalogue == 1)
        
        if ind == 1
            Start = []; ExpID_Start=[];
            UserVar2.type = "Inverse";
            UserVar2.home = pwd; 
            RunTablesToRead = uigetfile(".csv","Select run tables with experiments for startAGlen and startC","MultiSelect","on");
            for tt=1:numel(RunTablesToRead)
                UserVar2.Table = RunTablesToRead{tt};
                RunTable_Start = ANT_ReadWritetable(UserVar2,UserVar2.Table,[],'read');
                % scan table for ExpID, m and n                
                gaC_Start = RunTable_Start{:,'gaC'};
                gaA_Start = RunTable_Start{:,'gaA'};
                gsC_Start = RunTable_Start{:,'gsC'};
                gsA_Start = RunTable_Start{:,'gsA'};
                m_Start = RunTable_Start{:,'m'};
                n_Start = RunTable_Start{:,'n'};
                Finished = RunTable_Start{:,'Finished'};
                Start_tmp = [gaC_Start(:) gaA_Start(:) gsC_Start(:) gsA_Start(:) m_Start(:) n_Start(:)];
                Start_tmp = Start_tmp(Finished==1,:);
                Start = [Start; Start_tmp];
                ExpID_tmp = RunTable_Start{:,'ExpID'};
                ExpID_tmp = ExpID_tmp(Finished==1);
                ExpID_Start = [ExpID_Start; ExpID_tmp];            
            end   
        end
        [Idx,D] = knnsearch(Start,[gaC gaA gsC gsA m n],"Distance","seuclidean");
        startC = ExpID_Start(Idx);
        startAGlen = ExpID_Start(Idx);
        iterations = "5000+5000";
        spinupyears = "3";
        invertfor = "-logC-logA-";

    else % no start values
        startC = -9999;
        startAGlen = 0;
        iterations = "20";
        spinupyears = "0";
        invertfor = "-logC-";
    end

    Newrow = {'ANT_nsmbl',...               %Domain
        0,...                               %pgid
        0,...                               %ExpID
        0,...                               %Submitted
        "01/01/2000 00:00:00",...           %SubmissionTime
        0,...                               %Running
        0,...                               %Error
        "01/01/2000 00:00:00",...           %ErrorTime
        0,...                               %Finished
        "01/01/2000 00:00:00",...           %FinishedTime
        0,...                               %Restart
        iterations,...                      %InverseIterations
        0,...                               %InverseIterationsDone
        spinupyears,...                     %SpinupYears
        0,...                               %SpinupYearsDone
        invertfor,...                       %InvertFor
        GradientCalc,...                    %GradientCalc
        gsC,...                             %gsC
        gsA,...                             %gsA
        gaC,...                             %gaC
        gaA,...                             %gaA
        "-uv-dhdt-",...                     %Measurements
        dhdt_err,...                        %dhdt error scaling
        2000,...                            %Velocity
        2000,...                            %startGeometry
        startMesh,...                       %startMesh
        SlidingLaw,...                      %SlidingLaw
        m,...                               %m
        muk,...                             %muk
        priorC,...                          %priorC
        startC,...                          %startC
        n,...                               %n
        priorAGlen,...                      %priorAGlen
        startAGlen,...                      %startAGlen
        ""                                  %Comments
        };

    RunTable = [RunTable; Newrow];

end

[~] = ANT_ReadWritetable(UserVar,UserVar.Table,RunTable,'write');
for nn = 1:ceil(height(RunTable)/MaxTableLength)
    IndexRange = (nn-1)*MaxTableLength+1:min(nn*MaxTableLength,height(RunTable));
    RunTable_tmp = RunTable(IndexRange,:);
    TableName = erase(UserVar.Table,".csv")+"_"+string(nn)+".csv";
    [~] = ANT_ReadWritetable(UserVar,TableName,RunTable_tmp,'write');
end