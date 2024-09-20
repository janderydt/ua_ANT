function plot_PerturbationResults_Ensemble(diagnostic_to_plot,parameter_to_plot,cycles_to_plot)

addpath(getenv("froot_tools"));
addpath(getenv("froot_ua")+"cases/ANT");

if nargin==0
    diagnostic_to_plot = 'Delta_u'; % Delta_qGL, Delta_qOB, Delta_u
    parameter_to_plot = 'm'; % m, n, gaA, gaC, gsA, gsC
    cycles_to_plot = [1 2]; %[1 2]
    slidinglaws = ["Umbi"];
end

basins_to_analyze = {'F-G',...  % Getz
    'G-H',...  % PIG, Thwaites
    'H-Hp'}; % Abbott
file_with_perturbation_data_to_read = "perturbationdata_"+slidinglaws+".mat";
CtrlVar=Ua2D_DefaultParameters; 

%% load basins
filename = 'basins_IMBIE_v2.mat'; 
B = load(filename);
B = RemoveSmallIceRisesAndIslands(B);

%% prepare meshes
if diagnostic_to_plot=="Delta_u"

    UserVar.casefolder=pwd;
    UserVar.datafolder=pwd+"/../ANT_Data/";
    UserVar.Experiment="";
    BaseMesh='2000_2009_2014_2018_meshmin3000_meshmax100000_refined';
    UserVar = ANT_DefineBaseMesh(UserVar,BaseMesh);

    % 2000
    UserVar.Geometry=2000;
    UserVar=ANT_ApplyMeshModifications(UserVar);
    tmp = load(UserVar.InitialMeshFileName);
    MUA_2000 = tmp.MUA; 
    % identify basin id of each MUA node
    [MUA_2000.basins,~] = Define_Quantity_PerBasin(MUA_2000.coordinates(:,1),MUA_2000.coordinates(:,2),B,0);
    MUA_basinnames = erase({MUA_2000.basins(:).name},'-');
    basinnodes_all = [];
    for bb=1:numel(basins_to_analyze) 
        basin = char(erase(basins_to_analyze{bb},'-'));
        [~,BasinInd] = ismember(basin,MUA_basinnames);
        basinnodes = MUA_2000.basins(BasinInd).ind;
        basinnodes_all = [basinnodes_all; basinnodes];
    end
    ElementsToBeDeactivated=any(~ismember(MUA_2000.connectivity,basinnodes_all),2);
    [MUA_2000,MUA_2000.k,MUA_2000.l]=DeactivateMUAelements(CtrlVar,MUA_2000,ElementsToBeDeactivated);
    Ind_nan = find(isnan(MUA_2000.Boundary.x));
    if ~isempty(Ind_nan)
        MUA_2000.Boundary.x = MUA_2000.Boundary.x(1:Ind_nan(1)-1);
        MUA_2000.Boundary.y = MUA_2000.Boundary.y(1:Ind_nan(1)-1);
    end
    
    % 2018
    UserVar.Geometry=2018;
    UserVar=ANT_ApplyMeshModifications(UserVar);
    tmp = load(UserVar.InitialMeshFileName);
    MUA_2018 = tmp.MUA;
    [MUA_2018.basins,~] = Define_Quantity_PerBasin(MUA_2018.coordinates(:,1),MUA_2018.coordinates(:,2),B,0);
    MUA_basinnames = erase({MUA_2018.basins(:).name},'-'); 
    basinnodes_all=[];
    for bb=1:numel(basins_to_analyze)    
        basin = char(erase(basins_to_analyze{bb},'-'));
        [~,BasinInd] = ismember(basin,MUA_basinnames);
        basinnodes = MUA_2018.basins(BasinInd).ind;
        basinnodes_all = [basinnodes_all; basinnodes];
    end
    ElementsToBeDeactivated=any(~ismember(MUA_2018.connectivity,basinnodes_all),2);
    [MUA_2018,MUA_2018.k,MUA_2018.l]=DeactivateMUAelements(CtrlVar,MUA_2018,ElementsToBeDeactivated);
    Ind_nan = find(isnan(MUA_2018.Boundary.x));
    if ~isempty(Ind_nan)
        MUA_2018.Boundary.x = MUA_2018.Boundary.x(1:Ind_nan(1)-1);
        MUA_2018.Boundary.y = MUA_2018.Boundary.y(1:Ind_nan(1)-1);
    end

    delete(UserVar.InitialMeshFileName);
    delete(UserVar.MeshBoundaryCoordinatesFile);
end

%% load data
if exist(file_with_perturbation_data_to_read,"file")
    load(file_with_perturbation_data_to_read,"data","perturbation_experiments_analyzed","GF_2000","GF_2018");
    original_node_numbers = MUA_2000.k(find(~isnan(MUA_2000.k)));
    GF_2000.node=GF_2000.node(original_node_numbers);
    original_node_numbers = MUA_2018.k(find(~isnan(MUA_2018.k)));
    GF_2018.node=GF_2018.node(original_node_numbers);
else
    error(file_with_perturbation_data_to_read+" does not exist");
end
tmp = load("inversiondata_"+slidinglaws+".mat");
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
qGL_orig.total = zeros(numel(data),2); qOB_orig.total = zeros(numel(data),2);
qGL_calv.total = zeros(numel(data),2); qOB_calv.total = zeros(numel(data),2);
qGL_dhIS.total = zeros(numel(data),2); qOB_dhIS.total = zeros(numel(data),2);
qGL_dh.total = zeros(numel(data),2); qOB_dh.total = zeros(numel(data),2);
qGL_calvdh.total = zeros(numel(data),2); qOB_calvdh.total = zeros(numel(data),2);

misfit = zeros(numel(data),2);
gsA = zeros(numel(data),1);
gsC = zeros(numel(data),1);
gaA = zeros(numel(data),1);
gaC = zeros(numel(data),1);

Delta_u = []; Fu_original = [];

for ii=1:numel(data)

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
                qOB_orig.(basin)(ii,:) = data(ii).Original.qOB.(basin)(:)';
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

                for ff=["Calv","dhIS","dh","Calv_dh"]
                    
                    if contains(ff,"Calv")
                        MUA_target = MUA_2018;
                    else
                        MUA_target = MUA_2000;
                    end

                    for cc=1:size(data(ii).Original.speed,2)
                            
                        if contains(ff,"Calv")
                            original_node_numbers = MUA_2000.k(find(~isnan(MUA_2000.k)));
                            % interpolate original and perturbed speed to same grid
                            Ind_out = find(~inpoly2(MUA_target.coordinates,[MUA_2000.Boundary.x MUA_2000.Boundary.y]));
                            if isempty(Fu_original)
                                Fu_original = scatteredInterpolant(MUA_2000.coordinates(:,1),MUA_2000.coordinates(:,2),data(ii).Original.speed(original_node_numbers,cc),"natural");
                            else
                                Fu_original.Values = data(ii).Original.speed(original_node_numbers,cc);
                            end
                            u_original_interp = Fu_original(MUA_target.coordinates(:,1),MUA_target.coordinates(:,2));
                            u_original_interp(Ind_out) = nan;
                            original_node_numbers = MUA_target.k(find(~isnan(MUA_target.k)));
                            if size(data(ii).(char(ff)).speed,2)>=cc
                                du = data(ii).(char(ff)).speed(original_node_numbers,cc)-u_original_interp;
                            else
                                du = nan*original_node_numbers;
                            end
                        else
                            original_node_numbers = MUA_target.k(find(~isnan(MUA_target.k)));
                            if size(data(ii).(char(ff)).speed,2)>=cc
                                du = data(ii).(char(ff)).speed(original_node_numbers,cc)-data(ii).Original.speed(original_node_numbers,cc);
                            else
                                du = nan*original_node_numbers;
                            end
                        end
                        Intdu=FEintegrate2D(CtrlVar,MUA_target,du);
                        IntA=FEintegrate2D(CtrlVar,MUA_target,0*du+1);
                        Ind_notnan=find(~isnan(Intdu));
                        Delta_u.(char(ff)).(basin)(ii,cc) = sum(Intdu(Ind_notnan))/sum(IntA(Ind_notnan));
                        if bb==numel(basins_to_analyze)
                            if ~isfield(Delta_u.(char(ff)),'map')
                                Delta_u.(char(ff)).map = zeros(numel(data),MUA_target.Nnodes,size(data(ii).Original.speed,2));
                                %Delta_u.(char(ff)).mapcounter = 0;
                                %Delta_u.(char(ff)).mapnodes(:,cc) = basinnodes_all(:);
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
        gsA(ii,1) = [data(ii).Inverse.gsA];
        gsC(ii,1) = [data(ii).Inverse.gsC];
        gaA(ii,1) = [data(ii).Inverse.gaA];
        gaC(ii,1) = [data(ii).Inverse.gaC];

    end

    fprintf("Done %s out of %s.\n",string(ii),string(numel(data)));

end

m = [data(:).m];
n = [data(:).n];

if diagnostic_to_plot=="Delta_u"
    save("Delta_u.mat", "Delta_u","MUA_2000","MUA_2018","misfit","gsA","gsC","gaA","gaC","m","n");
end

CtrlVar.PlotXYscale = 1e3;

N = misfit-min(misfit);
N = N./max(N);
alphavalue = 1-N;
markersize = 50*alphavalue+eps;
newcolors = [0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54];

%% PLOTTING
figure; hold on;

switch diagnostic_to_plot

    case {'Delta_qGL','Delta_qOB'}

        H=fig('units','inches','width',70*12/72.27,'height',50*12/72.27,'fontsize',14,'font','Helvetica');

        ax1=subplot("position",[0.1 0.1 0.5 0.85]); hold on;
        ax2=subplot("position",[0.65 0.1 0.3 0.85]); hold on;

        switch parameter_to_plot
            case 'n'
                xdata = n(:);
            case 'm'
                xdata = m(:);
            case 'gsA'
                xdata = gsA(:);
            case 'gsC'
                xdata = gsC(:);
            case 'gaA'
                xdata = gaA(:);
            case 'gaC'
                xdata = gaC(:);
            otherwise
                xdata = [];
        end

        for cc=cycles_to_plot
        
            switch diagnostic_to_plot
        
                case 'Delta_qGL'
        
                    qorig = qGL_orig.total(:,cc);
                    ydata1 = qGL_calv.total(:,cc)-qorig;
                    ydata2 = qGL_dhIS.total(:,cc)-qorig;
                    ydata3 = qGL_dh.total(:,cc)-qorig;
                    ydata4 = qGL_calvdh.total(:,cc)-qorig;
                    yaxislabel = "\Deltaq_{GL} [Gt/yr]";
                    dq = Calc_Delta_qGL(cc,basins_to_analyze);
        
                case 'Delta_qOB'
        
                    qorig = qOB_orig.total(:,cc);
                    ydata1 = qOB_calv.total(:,cc)-qorig;
                    ydata2 = qOB_dhIS.total(:,cc)-qorig;
                    ydata3 = qOB_dh.total(:,cc)-qorig;
                    ydata4 = qOB_calvdh.total(:,cc)-qorig;
                    yaxislabel = "\Delta q_{OB} [Gt/yr]";
                    dq = nan;
        
                otherwise
        
                    continue;
        
            end
        
            if cc==1
                marker='o';                        
            elseif cc==2
                marker='s';               
            end

            scat((cc-1)*4+1)=scatter(ax1,xdata,ydata1,markersize(:,1),marker,"filled",'MarkerEdgeColor','none');
            scat((cc-1)*4+2)=scatter(ax1,xdata,ydata2,markersize(:,1),marker,"filled",'MarkerEdgeColor','none');
            scat((cc-1)*4+3)=scatter(ax1,xdata,ydata3,markersize(:,1),marker,"filled",'MarkerEdgeColor','none');
            scat((cc-1)*4+4)=scatter(ax1,xdata,ydata4,markersize(:,1),marker,"filled",'MarkerEdgeColor','none');

        end

        colororder(ax1,[newcolors;newcolors]);

        xlim(ax1,[floor(min(xdata)) ceil(max(xdata))]);
        
        for ii=0:4*cycles_to_plot(end)-1
            scat(ii+1).AlphaData = alphavalue(:,floor(ii/4)+1);
            scat(ii+1).MarkerFaceAlpha='flat';
        end
        
        for ii=1:4
            scat(ii)=plot(ax1,0,0,'o','color',newcolors(ii,:),'MarkerSize',10,'MarkerFaceColor',newcolors(ii,:));
        end

        if numel(cycles_to_plot)==2   
            scat(9)=plot(ax1,0,0,'ok','markersize',10);
            scat(10)=plot(ax1,0,0,'sk','markersize',10);
            legend(ax1,[scat(1:4) scat(9) scat(10)],{'Calving','Ice Shelf thickness','Ice thickness','Calving + Ice thickness','no spinup','with spinup'},...
            'NumColumns',2,'Location','northwest');
        else
            legend(ax1,scat(1:4),{'Calving','Ice Shelf thickness','Ice thickness','Calving + Ice thickness'},...
            'NumColumns',2,'Location','northwest');
        end

        xlabel(ax1,parameter_to_plot); ylabel(ax1,yaxislabel);
        
        if ismember(diagnostic_to_plot,["gsA","gsC","gaA","gaC"])
            ax1.XScale='log';
        end
        ax1.YScale='log';
        grid(ax1,"on")
        box(ax1,"on");
        ax1_ylim = ylim(ax1);
        yticks(ax1,[10:10:100 200:100:1000]);
        yticklabels(ax1,["10","","","","","","","","","100","","","","","","","","","1000"]);
        set(ax1,"YMinorGrid","off");
        title(ax1,"Ua perturbation experiments");

        % histogram
        [N,edges] = histcounts(dq,25,"Normalization","percentage");
        barh(ax2,0.5*(edges(1:end-1)+edges(2:end)),N);
        h=plot(ax2,[0 100],[100 100],'--k','linewidth',1.5); %change in dicharge (Gt/yr) according to Rignot 2019
        ax2.YScale='log';
        grid(ax2,"on")
        box(ax2,"on");
        xlabel(ax2,"Percentage");
        yticklabels(ax2,"");
        ylim(ax2,ax1_ylim);
        xlim(ax2,[0,15])
        yticks(ax2,[10:10:100 200:100:1000]);
        set(ax2,"YMinorGrid","off");
        legend(ax2,h,"Rignot et al. 2019","location","southeast");
        title(ax2,"'Measured' change");

        % Save
        pos = get(H,"Position");
        set(H,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);
        fname = "./Figures/"+diagnostic_to_plot+"_"+parameter_to_plot;
        print(H,fname,"-dpng","-r400");


    case "Delta_u"
        
        fields_to_plot = fields(Delta_u);

        % observed change in velocity
        load("../ANT_Data/ANT_Interpolants/GriddedInterpolants_1996-2003_MeaSUREs_ITSLIVE_Velocities_EXTRUDED","Fus","Fvs","Fxerr","Fyerr");
        ux2000 = Fus(MUA_2018.coordinates(:,1),MUA_2018.coordinates(:,2));
        uxerr2000 = Fxerr(MUA_2018.coordinates(:,1),MUA_2018.coordinates(:,2));
        uy2000 = Fvs(MUA_2018.coordinates(:,1),MUA_2018.coordinates(:,2));
        uyerr2000 = Fyerr(MUA_2018.coordinates(:,1),MUA_2018.coordinates(:,2));
        u2000 = hypot(ux2000,uy2000);
        uerr2000 = sqrt((ux2000.^2.*uxerr2000.^2+uy2000.^2.*uyerr2000.^2)./(ux2000.^2+uy2000.^2+eps));
        load("../ANT_Data/ANT_Interpolants/GriddedInterpolants_2018-2019_MeaSUREs_ITSLIVE_Velocities_EXTRUDED","Fus","Fvs","Fxerr","Fyerr");
        ux2018 = Fus(MUA_2018.coordinates(:,1),MUA_2018.coordinates(:,2));
        uxerr2018 = Fxerr(MUA_2018.coordinates(:,1),MUA_2018.coordinates(:,2));
        uy2018 = Fvs(MUA_2018.coordinates(:,1),MUA_2018.coordinates(:,2));
        uyerr2018 = Fyerr(MUA_2018.coordinates(:,1),MUA_2018.coordinates(:,2));
        u2018 = hypot(ux2018,uy2018);
        uerr2018 = sqrt((ux2018.^2.*uxerr2018.^2+uy2018.^2.*uyerr2018.^2)./(ux2018.^2+uy2018.^2+eps));

        dspeed_meas = u2018-u2000; %dspeed(dspeed<0)=nan;
        dspeed_meas_err = hypot(uerr2000,uerr2018);       

        for cc=cycles_to_plot

            %% plot distribution of misfit between modelled perturbation in speed and measured change in speed
            % prep data

            dspeed_meas_clean = dspeed_meas; 
            dspeed_meas_clean(abs(dspeed_meas)<1 | dspeed_meas<=0) = nan;

            mask = ones(MUA_2018.Nnodes,1);

            area_floating = FEintegrate2D(CtrlVar,MUA_2018,1-GF_2018.node); 
            area_grounded = FEintegrate2D(CtrlVar,MUA_2018,GF_2018.node); 

            for nn=1:size(Delta_u.Calv_dh.map,1)
                dspeed_mod = squeeze(Delta_u.Calv_dh.map(nn,:,cc))';
                dspeed_mod_clean = dspeed_mod;
                dspeed_mod_clean(abs(dspeed_mod)<1 | dspeed_mod<=0) = nan;

                dspeed_log_err = 1./(dspeed_mod_clean*log(10)).*dspeed_meas_err;
                
                misfit_tmp = abs(log10(dspeed_mod_clean) - log10(dspeed_meas_clean));%./(dspeed_log_err+eps)).^2;

                x = MUA_2018.coordinates(:,1);
                y = MUA_2018.coordinates(:,2);
                PIG_interestregion = find(x>-1.59e6 & x<-1.57e6 & y>-2.3e5 & y<-2.14e5);
                misfit_tmp_PIG(nn) = sum(misfit_tmp(PIG_interestregion));

                misfit_floating_tmp = FEintegrate2D(CtrlVar,MUA_2018,(1-GF_2018.node).*misfit_tmp);
                misfit_grounded_tmp = FEintegrate2D(CtrlVar,MUA_2018,GF_2018.node.*misfit_tmp);

                misfit_floating(nn) = sum(misfit_floating_tmp,"all","omitmissing")/sum(area_floating,"all","omitmissing");
                misfit_grounded(nn) = sum(misfit_grounded_tmp,"all","omitmissing")/sum(area_grounded,"all","omitmissing");

                disp("done "+num2str(nn)+" out of "+num2str(size(Delta_u.Calv_dh.map,1)));
            end
            
            H=fig('units','inches','width',120*12/72.27,'height',60*12/72.27,'fontsize',14,'font','Helvetica');


            figure; scatter([data(:).m],misfit_tmp_PIG)
            figure; scatter([data(:).m],misfit_grounded)
            figure; scatter([data(:).n],misfit_grounded)
            figure; scatter([data(:).m],misfit_floating)
            figure; scatter([data(:).n],misfit_floating)
            
            

            H=fig('units','inches','width',120*12/72.27,'height',60*12/72.27,'fontsize',14,'font','Helvetica');

            tlo(cc*999)=tiledlayout(H,2,numel(fields_to_plot)+1,'TileSpacing','none','TileIndexing', 'columnmajor');

            for ff=1:numel(fields_to_plot)

                deltau_av = mean(Delta_u.(fields_to_plot{ff}).map(:,:,cc),1,"omitmissing");
                deltau_std = std(Delta_u.(fields_to_plot{ff}).map(:,:,cc),1,"omitmissing");

                if contains(fields_to_plot{ff},'Calv')
                    MUA = MUA_2018;
                    MUA2 = MUA_2000;
                else
                    MUA = MUA_2000;
                    MUA2 = MUA_2018;
                end
                xmin = min(MUA.coordinates(:,1)); xmax = max(MUA.coordinates(:,1));
                ymin = min(MUA.coordinates(:,2)); ymax = max(MUA.coordinates(:,2));

                %% ensemble-average change in speed
                ax(ff)=nexttile(tlo(cc*999)); hold on;
                log_deltau_av = deltau_av(:);
                log_deltau_av(log_deltau_av<0)=nan;
                log_deltau_av = log10(log_deltau_av);
                %log_deltau_av(isnan(log_deltau_av))
                PlotNodalBasedQuantities_JDR(ax(ff),MUA.connectivity,MUA.coordinates,log_deltau_av(:),CtrlVar);
                plot(MUA.Boundary.x/CtrlVar.PlotXYscale,MUA.Boundary.y/CtrlVar.PlotXYscale,'-k');
                plot(MUA2.Boundary.x/CtrlVar.PlotXYscale,MUA2.Boundary.y/CtrlVar.PlotXYscale,'--k');
                PlotGroundingLines(CtrlVar,MUA_2018,GF_2018,[],[],[],'-k','linewidth',1);
                CM = flipdim(othercolor('RdYlBu8',15),1);
                %CM1 = othercolor('RdYlBu8',3);
                %CM = CM([15 19:40],:);
                colormap(ax(ff),CM);
                title(ax(ff),fields_to_plot{ff},"Interpreter","none");
                axis(ax(ff),"off");
                caxis(ax(ff),[-0.25 3.5]);
                xlim(ax(ff),[xmin xmax]/CtrlVar.PlotXYscale);
                ylim(ax(ff),[ymin ymax]/CtrlVar.PlotXYscale);
                axis(ax(ff),"equal");

                %% standard deviation 
                ax(numel(fields_to_plot)+ff)=nexttile(tlo(cc*999)); hold on;
                PlotNodalBasedQuantities_JDR(ax(numel(fields_to_plot)+ff),MUA.connectivity,MUA.coordinates,log10(deltau_std(:)+eps),CtrlVar);
                plot(MUA.Boundary.x/CtrlVar.PlotXYscale,MUA.Boundary.y/CtrlVar.PlotXYscale,'-k');
                plot(MUA2.Boundary.x/CtrlVar.PlotXYscale,MUA2.Boundary.y/CtrlVar.PlotXYscale,'--k');
                PlotGroundingLines(CtrlVar,MUA_2018,GF_2018,[],[],[],'-k','linewidth',1);
                %CM = turbo(13); CM(1,:)=[1 1 1];
                colormap(ax(numel(fields_to_plot)+ff),CM);
                axis(ax(numel(fields_to_plot)+ff),"off");
                caxis(ax(numel(fields_to_plot)+ff),[-0.25 3.5]);
                xlim(ax(numel(fields_to_plot)+ff),[xmin xmax]/CtrlVar.PlotXYscale);
                ylim(ax(numel(fields_to_plot)+ff),[ymin ymax]/CtrlVar.PlotXYscale);
                axis(ax(numel(fields_to_plot)+ff),"equal");

                if ff==numel(fields_to_plot)
                     cb=colorbar(ax(2*numel(fields_to_plot))); cb.Label.String="Ensemble standard deviation of log_{10}(\Deltau) [m/yr]";
                end

            end          

            ax(numel(fields_to_plot)+1)=nexttile(tlo(cc*999)); hold on;
            PlotNodalBasedQuantities_JDR(ax(numel(fields_to_plot)+1),MUA_2018.connectivity,MUA_2018.coordinates,log10(dspeed_meas+eps),CtrlVar);
            plot(MUA_2018.Boundary.x/CtrlVar.PlotXYscale,MUA_2018.Boundary.y/CtrlVar.PlotXYscale,'-k');
            plot(MUA_2000.Boundary.x/CtrlVar.PlotXYscale,MUA_2000.Boundary.y/CtrlVar.PlotXYscale,'--k');
            PlotGroundingLines(CtrlVar,MUA_2018,GF_2018,[],[],[],'-k','linewidth',1);
            %CM = turbo(13); CM(1,:)=[1 1 1];
            colormap(ax(numel(fields_to_plot)+1),CM);
            title(ax(numel(fields_to_plot)+1),"Observations");
            axis(ax(numel(fields_to_plot)+1),"off");
            caxis(ax(numel(fields_to_plot)+1),[-0.25 3.5]);
            xlim(ax(numel(fields_to_plot)+1),[xmin xmax]/CtrlVar.PlotXYscale);
            ylim(ax(numel(fields_to_plot)+1),[ymin ymax]/CtrlVar.PlotXYscale);
            axis(ax(numel(fields_to_plot)+1),"equal");

            cb2=colorbar(ax(numel(fields_to_plot)+1)); cb2.Label.String="Log_{10} change in speed [m/yr]";

            title(tlo(cc*999),"Cycle "+num2str(cc));

            % Save
            pos = get(H,"Position");
            set(H,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);
            fname = "./Figures/"+diagnostic_to_plot+"cycle"+string(cc);
            print(H,fname,"-dpng","-r400");


        end

    otherwise 

        error("unknown case")
       
end


