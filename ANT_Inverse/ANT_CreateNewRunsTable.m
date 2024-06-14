function ANT_CreateNewRunsTable(X,GradientCalc,SlidingLaw)

addpath('../');

UserVar.home = pwd;
UserVar.Table = "NewRuns_"+GradientCalc+".csv";
if ~exist(UserVar.Table,"file")
    copyfile("EmptyTable.csv",UserVar.Table);
end
UserVar.type = "Inverse";

RunTable = ANT_ReadWritetable(UserVar,UserVar.Table,[],'read');

for ind=1:size(X,1)

    startMesh = "2000_2009_2014_2018_meshmin3000_meshmax100000_refined";
    
    taub = 80; m = round(X.m(ind)*100)/100; ub = 100;%X(ind,6);
    priorC = ub/taub^m;
    muk = 0.5;
    
    tau = 80; n = round(X.n(ind)*100)/100; eps = 0.0026;%X(ind,8);
    priorAGlen = eps/tau^n;
    
    % If GradientCalc=Adjoint, find Fixpoint inversion with nearest (m,n)
    % pair as start value. The code will rescale C and AGlen depending on m and n 
    if GradientCalc == "Adjoint"
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
        iterations,...                            %InverseIterations
        0,...                               %InverseIterationsDone
        spinupyears,...                             %SpinupYears
        0,...                               %SpinupYearsDone
        invertfor,...                        %InvertFor
        GradientCalc,...                    %GradientCalc
        round(X.gsC(ind)*10)/10,...         %gsC
        round(X.gsA(ind)*10)/10,...         %gsA
        round(X.gaC(ind)*10)/10,...         %gaC
        round(X.gaA(ind)*10)/10,...         %gaA
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
RunTable
[~] = ANT_ReadWritetable(UserVar,UserVar.Table,RunTable,'write');
