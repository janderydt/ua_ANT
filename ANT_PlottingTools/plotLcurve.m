function plotLcurve

ExpID = [1814 1278 1546 1792 1959 1655 1035 1849];

cm = crameri('roma',numel(ExpID));

for i = 1:numel(ExpID)

    load("../ANT_Inverse/ANT_Inverse_"+string(ExpID(i))+"/ANT_Inverse_"+string(ExpID(i))+"-RestartFile.mat");

    if contains(CtrlVarInRestartFile.Inverse.Regularize.Field,"log")
        yaC(i) = CtrlVarInRestartFile.Inverse.Regularize.logC.ga;
        yaA(i) = CtrlVarInRestartFile.Inverse.Regularize.logAGlen.ga;
        ysC(i) = CtrlVarInRestartFile.Inverse.Regularize.logC.gs;
        ysA(i) = CtrlVarInRestartFile.Inverse.Regularize.logAGlen.gs;
    else
        yaC(i) = CtrlVarInRestartFile.Inverse.Regularize.C.ga;
        yaA(i) = CtrlVarInRestartFile.Inverse.Regularize.AGlen.ga;
        ysC(i) = CtrlVarInRestartFile.Inverse.Regularize.C.gs;
        ysA(i) = CtrlVarInRestartFile.Inverse.Regularize.AGlen.gs;
    end

    I(i) = InvFinalValues.I;
    R(i) = InvFinalValues.R/(ysC(i)^2);

    figure(111); hold on; 
    %yyaxis left
    g(i)=plot(RunInfo.Inverse.Iterations,RunInfo.Inverse.J,'-x','LineWidth',2,'color',cm(i,:));
    ylabel('J','interpreter','latex');
    
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
    xlabel('Inverse iteration','interpreter','latex');
    hold off

end

figure(111); legend(g(:),string(ysC))

figure; semilogx(R,I,'ok');
for i=1:numel(ExpID)
    text(R(i),I(i),num2str(ysC(i)))
end


