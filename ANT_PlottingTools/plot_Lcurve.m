function plot_Lcurve

UserVar.cycle = 1;

%ExpID =[1849 1546 1035 1655 1959 1814 1792 1278];
RunTable = "RunTable_ARCHER2_"+string([2 5 8])+".csv";

gaA=[]; gaC=[]; gsA=[]; gsC=[]; m=[]; n=[]; finished=[]; error=[];
Ind_finished=[]; Ind_error=[]; ExpID=[];

for tt=1:numel(RunTable)

    UserVar.home = pwd; 
    UserVar.type = "Inverse";
    UserVar.Table = UserVar.home+"/../ANT_"+UserVar.type+"/"+RunTable(tt);

    RunTable_tmp = ANT_ReadWritetable(UserVar,UserVar.Table,[],'read');

    gaA = [gaA; RunTable_tmp{:,'gaA'}];
    gaC = [gaC; RunTable_tmp{:,'gaC'}];
    gsA = [gsA; RunTable_tmp{:,'gsA'}];
    gsC = [gsC; RunTable_tmp{:,'gsC'}];
    m = [m; RunTable_tmp{:,'m'}];
    n = [n; RunTable_tmp{:,'n'}];
    finished_tmp = RunTable_tmp{:,'Finished'};
    finished = [finished; finished_tmp];
    error_tmp = RunTable_tmp{:,'Error'};
    error = [error; error_tmp];
    ExpID_tmp = RunTable_tmp{:,'ExpID'};
    ExpID = [ExpID; ExpID_tmp];

    Ind_finished_tmp = find(finished_tmp==1);
    Ind_finished = [Ind_finished; Ind_finished_tmp];
    Ind_error_tmp = find(error_tmp==1);
    Ind_error = [Ind_error; Ind_error_tmp];

end

gaA = gaA(Ind_finished);
gaC = gaC(Ind_finished);
gsA = gsA(Ind_finished);
gsC = gsC(Ind_finished);
m = m(Ind_finished);
n = n(Ind_finished);
ExpID = ExpID(Ind_finished);

cm = crameri('roma',numel(ExpID));

for i = 1:numel(ExpID)

    load("../ANT_Inverse/cases/ANT_nsmbl_Inverse_"+string(ExpID(i))+"/ANT_nsmbl_Inverse_"+string(ExpID(i))+...
        "-RestartFile_InverseCycle"+string(UserVar.cycle)+".mat",...
        "InvFinalValues");

    I(i) = InvFinalValues.I;
    R_gsC(i) = InvFinalValues.R/(gsC(i)^2);
    R_gsA(i) = InvFinalValues.R/(gsA(i)^2);
    R_gaC(i) = InvFinalValues.R/(gaC(i)^2);
    R_gaA(i) = InvFinalValues.R/(gaA(i)^2);

    %figure(111); hold on; 
    %yyaxis left
    %g(i)=plot(RunInfo.Inverse.Iterations,RunInfo.Inverse.J,'-x','LineWidth',2,'color',cm(i,:));
    %ylabel('J','interpreter','latex');
    
%     if ~isempty(RunInfo.Inverse.GradNorm)  && ~all(isnan(RunInfo.Inverse.GradNorm)) ...
%             &&  numel(RunInfo.Inverse.Iterations) == numel(RunInfo.Inverse.GradNorm)
% 
%         hold off
%         yyaxis right
%         semilogy(RunInfo.Inverse.Iterations,RunInfo.Inverse.GradNorm,'-r+')
%         ylabel('Norm of gradient ($l_2$)','interpreter','latex')
%         legend('Objective function','$\| \hbox{Gradient} \|$','Location','northeast','interpreter','latex')
%         
%     end
    
%     yyaxis left
    %xlabel('Inverse iteration','interpreter','latex');
    %hold off

    disp("Done "+string(i)+" out of "+string(numel(ExpID)));

end

%figure(111); legend(g(:),string(ysC))

figure; scatter(log10(R_gsC),I,20,gsC,'filled'); title("\gamma_sC"); grid on; box on; xlabel("log_{10}(R/\gamma_{sC}^2)"); ylabel("I");
figure; scatter(log10(R_gsA),I,20,gsA,'filled'); title("\gamma_sA"); grid on; box on; xlabel("log_{10}(R/\gamma_{sA}^2)"); ylabel("I");
figure; scatter(log10(R_gaC),I,20,gaC,'filled'); title("\gamma_aC"); grid on; box on; xlabel("log_{10}(R/\gamma_{aC}^2)"); ylabel("I");
figure; scatter(log10(R_gaA),I,20,gaA,'filled'); title("\gamma_aA"); grid on; box on; xlabel("log_{10}(R/\gamma_{aA}^2)"); ylabel("I");
% for i=1:numel(ExpID)
%     text(R(i),I(i),num2str(ysC(i)))
% end


