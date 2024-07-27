function plot_PerturbationResults_Ensemble

addpath(getenv("froot_tools"));

diagnostic_to_plot = 'Delta_u'; % Delta_qGL, Delta_qOB, Delta_u
basins_to_analyze = {'F-G',...  % Getz
    'G-H',...  % PIG, Thwaites
    'H-Hp'}; % Abbot 
file_with_perturbation_data_to_read = "perturbationdata_intermediate.mat";

%% load data
if exist(file_with_perturbation_data_to_read,"file")
    load(file_with_perturbation_data_to_read);
else
    error(file_with_perturbation_data_to_read+" does not exist");
end
tmp = load("inversiondata.mat");
data_inverse = tmp.data;

% available drainage basins
available_basins = fieldnames(data(1).Original.qGL);
% check that data for basins_to_analyze is available
for bb=1:numel(basins_to_analyze)
    if ~ismember(erase(basins_to_analyze{bb},'-'),available_basins)
        error("Data for basin "+basins_to_analyze{bb}+" is not available.")
    end
end

%% gather data in a user-friendly format
qGL_orig.total = zeros(numel(data),10); qOB_orig.total = zeros(numel(data),10);
qGL_calv.total = zeros(numel(data),10); qOB_calv.total = zeros(numel(data),10);
qGL_dhIS.total = zeros(numel(data),10); qOB_dhIS.total = zeros(numel(data),10);
qGL_dh.total = zeros(numel(data),10); qOB_dh.total = zeros(numel(data),10);
qGL_calvdh.total = zeros(numel(data),10); qOB_calvdh.total = zeros(numel(data),10);
Delta_u = []; CtrlVar=Ua2D_DefaultParameters; ElementsToBeDeactivated=[];

for ii=1:numel(data)-1

    basinnodes_all = [];

    for bb=1:numel(basins_to_analyze)
    
        basin = char(erase(basins_to_analyze{bb},'-'));

        switch diagnostic_to_plot

            % grounding line flux
            case "Delta_qGL"
                qGL_orig.(basin)(ii,:) = data(ii).Original.qGL.(basin)(:)';
                qGL_calv.(basin)(ii,:) = data(ii).Calv.qGL.(basin)(:)';
                qGL_dhIS.(basin)(ii,:) = data(ii).dhIS.qGL.(basin)(:)';
                qGL_dh.(basin)(ii,:) = data(ii).dh.qGL.(basin)(:)';
                qGL_calvdh.(basin)(ii,:) = data(ii).Calv_dh.qGL.(basin)(:)';
                qGL_orig.total(ii,:) = qGL_orig.total(ii,:)+qGL_orig.(basin)(ii,:);
                qGL_calv.total(ii,:) = qGL_calv.total(ii,:)+qGL_calv.(basin)(ii,:);
                qGL_dhIS.total(ii,:) = qGL_dhIS.total(ii,:)+qGL_dhIS.(basin)(ii,:);
                qGL_dh.total(ii,:) = qGL_dh.total(ii,:)+qGL_dh.(basin)(ii,:);
                qGL_calvdh.total(ii,:) = qGL_calvdh.total(ii,:)+qGL_calvdh.(basin)(ii,:);

            % open boundary (calving) flux    
            case "Delta_qOB"
                qOB.orig.(basin)(ii,:) = data(ii).Original.qOB.(basin)(:)';
                qOB_calv.(basin)(ii,:) = data(ii).Calv.qOB.(basin)(:)';
                qOB_dhIS.(basin)(ii,:) = data(ii).dhIS.qOB.(basin)(:)';
                qOB_dh.(basin)(ii,:) = data(ii).dh.qOB.(basin)(:)';
                qOB_calvdh.(basin)(ii,:) = data(ii).Calv_dh.qOB.(basin)(:)';
                qOB_orig.total(ii,:) = qOB_orig.total(ii,:)+qOB_orig.(basin)(ii,:);
                qOB_calv.total(ii,:) = qOB_calv.total(ii,:)+qOB_calv.(basin)(ii,:);
                qOB_dhIS.total(ii,:) = qOB_dhIS.total(ii,:)+qOB_dhIS.(basin)(ii,:);
                qOB_dh.total(ii,:) = qOB_dh.total(ii,:)+qOB_dh.(basin)(ii,:);
                qOB_calvdh.total(ii,:) = qOB_calvdh.total(ii,:)+qOB_calvdh.(basin)(ii,:);

            
            % changes in speed - only keep data for the selected basins
            case "Delta_u"

                MUA_basinnames = erase({MUA_coarse.basins(:).name},'-');
                [~,BasinInd] = ismember(basin,MUA_basinnames);
                basinnodes = MUA_coarse.basins(BasinInd).ind;
                basinnodes_all = [basinnodes_all; basinnodes];

                for cc=1:size(data(ii).Original.speed,2)
                    for ff=["Calv","dhIS","dh","Calv_dh"]
                        du = data(ii).(char(ff)).speed(:,cc)-data(ii).Original.speed(:,cc);                        
                        tmp = NaN*du; tmp(basinnodes) = du(basinnodes);
                        du = tmp;
                        Intdu=FEintegrate2D([],MUA_coarse,du);
                        IntA=FEintegrate2D([],MUA_coarse,0*du+1);
                        Ind_notnan=find(~isnan(Intdu));
                        Delta_u.(char(ff)).(basin)(ii,cc) = sum(Intdu(Ind_notnan))/sum(IntA(Ind_notnan));
                        if bb==numel(basins_to_analyze)
                            du = data(ii).(char(ff)).speed(:,cc)-data(ii).Original.speed(:,cc);
                            if isempty(ElementsToBeDeactivated)
                                ElementsToBeDeactivated=any(~ismember(MUA_coarse.connectivity,basinnodes_all),2);
                                [MUA_new,k,~]=DeactivateMUAelements(CtrlVar,MUA_coarse,ElementsToBeDeactivated);
                            end
                            du = du(k);
                            Intdu=FEintegrate2D([],MUA_new,du);
                            IntA=FEintegrate2D([],MUA_new,0*du+1);
                            Ind_notnan=find(~isnan(Intdu));
                            Delta_u.(char(ff)).total(ii,cc) = sum(Intdu(Ind_notnan))/sum(IntA(Ind_notnan));
                            if ~isfield(Delta_u.(char(ff)),'map')
                                Delta_u.(char(ff)).map = zeros(numel(data),MUA_new.Nnodes,size(data(ii).Original.speed,2));
                                %Delta_u.(char(ff)).mapcounter = 0;
                                Delta_u.(char(ff)).mapnodes(:,cc) = basinnodes_all(:);
                            end
                            Delta_u.(char(ff)).map(ii,:,cc) = du(:);
                            %Delta_u.(char(ff)).mapcounter = Delta_u.(char(ff)).mapcounter+1;                    
                            %Delta_u.(char(ff)).map_max(:,cc)
                            %Delta_u.(char(ff)).map_min(:,cc)
                        end
                    end
                end
     
            otherwise

                error("diagnostic_to_plot "+diagnostic_to_plot+" not known.");

        end

        % inversion parameters
        misfit(ii,:)=data(ii).Inverse.misfit(:)';
        gsA(ii) = [data(ii).Inverse.gsA];
        gsC(ii) = [data(ii).Inverse.gsC];
        gaA(ii) = [data(ii).Inverse.gaA];
        gaC(ii) = [data(ii).Inverse.gaC];

    end

    fprintf("Done %s out of %s.\n",string(ii),string(numel(data)));

end

CtrlVar=Ua2D_DefaultParameters;
CtrlVar.PlotXYscale = 1e3;

A = gpuArray(zeros(size(Delta_u.Calv_dh.map(:,:,1))));
A = Delta_u.Calv_dh.map(:,:,1); A = A'; %rows: nodes, columns: experiments

nx = size(A,1); ny = size(A,2);

tic
[U,S,V] = svd(A,'econ');
toc

energy = cumsum(diag(S))/sum(diag(S));

figure(111), tlo1=tiledlayout(2,2,'TileSpacing','tight'); title(tlo1,'A');
nexttile;
imagesc(A), axis off; colormap(slanCM('YlGnBu')); cb1=colorbar;
title('Original');

figure(222), tlo2=tiledlayout(2,2,'TileSpacing','tight'); title(tlo2,'t = 1 year');
nexttile; 
PlotNodalBasedQuantities_JDR(gca,MUA_new.connectivity,MUA_new.coordinates,A(:,1),CtrlVar), axis equal, axis off, colormap(slanCM('YlGnBu')); cb2=colorbar(gca);
title('Original');

figure(333), tlo3=tiledlayout(2,2,'TileSpacing','tight'); title(tlo3,'t = 850 years');
nexttile; 
PlotNodalBasedQuantities_JDR(gca,MUA_new.connectivity,MUA_new.coordinates,A(:,end),CtrlVar), axis equal, axis off, colormap(slanCM('YlGnBu')); cb3=colorbar(gca);
title('Original');

figure(444); tlo4=tiledlayout(2,5,'TileSpacing','tight');
for ii=1:10
    nexttile;
    PlotNodalBasedQuantities_JDR(gca,MUA_new.connectivity,MUA_new.coordinates,U(:,ii),CtrlVar);
    colormap(othercolor('RdYlBu8'));
    title(['mode',num2str(ii)]); caxis([-0.02 0.02]);
    axis tight; axis off;
    %cb=colormap; cb.visible='off';
end
cb4=colorbar(gca);
plotind = 2;
for r = [10 50 100]
    Xapprox = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
    %Xapprox(Xapprox<0)=0; 
    %Xapprox(Xapprox>1)=1;
    figure(111); nexttile, imagesc(Xapprox), axis off, colormap(slanCM('YlGnBu'));
    figure(222); nexttile, PlotNodalBasedQuantities_JDR(gca,MUA_new.connectivity,MUA_new.coordinates,Xapprox(:,1),CtrlVar), axis equal, axis off, colormap(slanCM('YlGnBu'));
    title(['r=',num2str(r,'%d'),', ',num2str(100*energy(r),'%2.2f'),'% cumulative energy']);
    figure(333); nexttile, PlotNodalBasedQuantities_JDR(gca,MUA_new.connectivity,MUA_new.coordinates,Xapprox(:,end),CtrlVar), axis equal, axis off, colormap(slanCM('YlGnBu'));
    title(['r=',num2str(r,'%d'),', ',num2str(100*energy(r),'%2.2f'),'% cumulative energy']);
    plotind = plotind + 1;
end

cb1.Layout.Tile = 'east';
cb2.Layout.Tile = 'east';
cb3.Layout.Tile = 'east';
cb4.Layout.Tile = 'east';

%% singular values
figure, subplot(1,2,1)
semilogy(diag(S),'k','linewidth',2), grid on;
xlabel('r')
ylabel('singlar value, \sigma_r')
set(gca,'fontsize',14)
subplot(1,2,2)
plot(energy,'k','linewidth',2), grid on;
xlabel('r')
ylabel('cumulative energy')
set(gca,'fontsize',14)



Delta_u_mean=NaN*MUA_coarse.coordinates(:,1);  
Delta_u_std=NaN*MUA_coarse.coordinates(:,1);  
Delta_u_mean(Delta_u.Calv_dh.mapnodes(:,1))=mean(Delta_u.Calv_dh.map(:,:,1),1,"omitmissing");
Delta_u_std(Delta_u.Calv_dh.mapnodes(:,1))=std(Delta_u.Calv_dh.map(:,:,1),0,1,"omitmissing");

figure; PlotMeshScalarVariable([],MUA_coarse,Delta_u_mean);
figure; PlotMeshScalarVariable([],MUA_coarse,Delta_u_std);


m = [data(:).m];
n = [data(:).n];

N = misfit-min(misfit);
N = N./max(N);
alphavalue = 1-N;
markersize = 50*alphavalue+eps;
newcolors = [0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54];

switch variable_to_plot
    case 'n'
        xdata = n;
    case 'm'
        xdata = m;
    case 'gsA'
        xdata = gsA;
    case 'gsC'
        xdata = gsC;
    case 'gaA'
        xdata = gaA;
    case 'gaC'
        xdata = gaC;
end

%% Plotting
figure; hold on;

s(1)=scatter(xdata,calv(:,1)-orig(:,1),markersize(:,1),'o',"filled",'MarkerEdgeColor','none');
s(2)=scatter(xdata,dhIS(:,1)-orig(:,1),markersize(:,1),'o',"filled",'MarkerEdgeColor','none');
s(3)=scatter(xdata,dh(:,1)-orig(:,1),markersize(:,1),'o',"filled",'MarkerEdgeColor','none');
s(4)=scatter(xdata,calvdh(:,1)-orig(:,1),markersize(:,1),'o',"filled",'MarkerEdgeColor','none');

s(5)=scatter(xdata,calv(:,2)-orig(:,2),markersize(:,2),'s',"filled",'MarkerEdgeColor','none');
s(6)=scatter(xdata,dhIS(:,2)-orig(:,2),markersize(:,2),'s',"filled",'MarkerEdgeColor','none');
s(7)=scatter(xdata,dh(:,2)-orig(:,2),markersize(:,2),'s',"filled",'MarkerEdgeColor','none');
s(8)=scatter(xdata,calvdh(:,2)-orig(:,2),markersize(:,2),'s',"filled",'MarkerEdgeColor','none');

colororder([newcolors;newcolors]);

xlim([floor(min(xdata)) ceil(max(xdata))]);

for ii=0:7
    s(ii+1).AlphaData = alphavalue(:,floor(ii/4)+1);
    s(ii+1).MarkerFaceAlpha='flat';
end

for ii=1:4
    s(ii)=plot(0,0,'o','color',newcolors(ii,:),'MarkerSize',10,'MarkerFaceColor',newcolors(ii,:));
end
s(5)=plot(0,0,'ok','markersize',10);
s(6)=plot(0,0,'sk','markersize',10);

legend(s(1:6),{'Calving','Ice Shelf thickness','Ice thickness','Calving + Ice thickness','no spinup','with spinup'},...
    'NumColumns',2,'Location','northwest');

xlabel(variable_to_plot); ylabel('\Delta q_{GL} [Gt/yr]');
ax=gca;
ax.YScale='log';
ax.XScale='log';

grid on;
box on;


