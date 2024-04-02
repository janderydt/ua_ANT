function ANT_CreateNewRunsTable(X)

addpath('../');

UserVar.Table = 'NewRuns.csv';
UserVar.type = 'Inverse';

%% Fixpoint calculations, Unrefined mesh, No dh/dt, Weertman sliding
RunTable_fixpoint = ANT_ReadWritetable(UserVar,[],'read');

for ind=1:size(X,1)

    tau = 80; m = round(X(ind,5)*100)/100; ub = X(ind,6);
    priorC = ub/tau^m;
    n = round(X(ind,7)*100)/100; eps = X(ind,8);
    priorAGlen = eps/tau^n;

    Newrow = {0,...                         %pgid
        0,...                               %ExpID
        0,...                               %Submitted
        "01/01/2000 00:00:00",...           %SubmissionTime
        0,...                               %Running
        0,...                               %Error
        "01/01/2000 00:00:00",...           %ErrorTime
        0,...                               %Finished
        "01/01/2000 00:00:00",...           %FinishedTime
        0,...                               %Restart
        "10000+1000",...                    %InverseIterations
        0,...                               %InverseIterationsDone
        "1",...                             %SpinupYears
        0,...                               %SpinupYearsDone
        "-logC-",...                        %InvertFor
        "FixPoint",...                      %GradientCalc
        round(X(ind,1)*10)/10,...           %gsC
        round(X(ind,2)*10)/10,...           %gsA
        round(X(ind,3)*10)/10,...           %gaC
        round(X(ind,4)*10)/10,...           %gaA
        "-uv-",...                          %Measurements
        2000,...                            %Velocity
        2000,...                            %startGeometry
        "2000_2009_2014_2018_meshmin3000_meshmax100000",...      %startMesh
        "Weertman",...                      %SlidingLaw
        m,...                               %m
        priorC,...                          priorC
        0,...                               %startC
        n,...                               %n
        priorAGlen,...                      %priorAGlen
        0,...                               %startAGlen
        ""                                  %Comments
        };

    RunTable_fixpoint = [RunTable_fixpoint; Newrow];

end

[~] = ANT_ReadWritetable(UserVar,RunTable_fixpoint,'write');
