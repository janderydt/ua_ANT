function ANT_GenerateRunTable(X,GradientCalc,SlidingLaw,Enrich)

%% INPUTS:
% -> X is the runtable created by ANT_GenerateRunTable or
%    ANT_ExpandExistingRunTable
% -> GradientCalc is either 'fixpoint' or 'adjoint'
% -> SlidingLaw is any valid type of sliding law in Ua, e.g. Weertman
% -> Enrich is 0 or 1, it only affects the name of the csv file that is
%    produced

addpath('../');

UserVar.home = pwd;
UserVar.Table = "NewTable_"+GradientCalc+"_"+SlidingLaw+"_Enriched"+string(Enrich)+".csv";
UserVar.type = "Inverse";

if ~exist(UserVar.Table,"file")
    copyfile("EmptyTable.csv",UserVar.Table);
end

RunTable = ANT_ReadWritetable(UserVar,UserVar.Table,[],'read');

for ind=1:size(X,1)

    startMesh = "2000_2009_2014_2018_meshmin3000_meshmax100000_refined";
    
    taub = 80; m = round(X.m(ind)*100)/100; ub = 100;%X(ind,6);
    priorC = ub/taub^m;
    muk = 0.5;
    
    tau = 80; n = round(X.n(ind)*100)/100; eps = 0.0026;%X(ind,8);
    priorAGlen = eps/tau^n;

    gsC = round(X.gsC(ind)*10)/10;         %gsC
    gsA = round(X.gsA(ind)*10)/10;         %gsA
    gaC = round(X.gaC(ind)*10)/10;         %gaC
    gaA = round(X.gaA(ind)*10)/10;         %gaA
    
    % If GradientCalc=Adjoint, find Fixpoint inversion with nearest (m,n)
    % pair as start value. The code will rescale C and AGlen depending on m and n 
    if GradientCalc == "Adjoint" && Enrich == 0
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
        iterations = "10000+1000";
        spinupyears = "1";
        invertfor = "-logC-logA-";
    elseif GradientCalc == "Adjoint" && Enrich == 1
        if ind == 1
            UserVar2.type = "Inverse";
            UserVar2.home = pwd; 
            UserVar2.Table = input("Path of run table with experiments for startAGlen and startC: ","s");
            RunTable_Start = ANT_ReadWritetable(UserVar2,UserVar2.Table,[],'read');
            % scan table for ExpID, m and n
            ExpID_Start = RunTable_Start{:,'ExpID'};
            gaC_Start = RunTable_Start{:,'gaC'};
            gaA_Start = RunTable_Start{:,'gaA'};
            gsC_Start = RunTable_Start{:,'gsC'};
            gsA_Start = RunTable_Start{:,'gsA'};
            m_Start = RunTable_Start{:,'m'};
            n_Start = RunTable_Start{:,'n'};
            Start = [gaC_Start(:) gaA_Start(:) gsC_Start(:) gsA_Start(:) m_Start(:) n_Start(:)];
        end
        [Idx,D] = knnsearch(Start,[gaC gaA gsC gsA m n],"Distance","seuclidean");
        startC = ExpID_Start(Idx);
        startAGlen = ExpID_Start(Idx);
        iterations = "5000+1000";
        spinupyears = "1";
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
        "-uv-",...                          %Measurements
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
