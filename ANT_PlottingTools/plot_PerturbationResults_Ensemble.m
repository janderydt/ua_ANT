function plot_PerturbationResults_Ensemble(diagnostic_to_plot,parameter_to_plot,cycles_to_plot)

addpath(getenv("froot_tools"));

if nargin==0
    diagnostic_to_plot = 'Delta_u'; % Delta_qGL, Delta_qOB, Delta_u
    parameter_to_plot = 'm'; % m, n, gaA, gaC, gsA, gsC
    cycles_to_plot = [1 2]; %[1 2]
end

basins_to_analyze = {'F-G',...  % Getz
    'G-H',...  % PIG, Thwaites
    'H-Hp'}; % Abbot 
file_with_perturbation_data_to_read = "perturbationdata_tmp.mat";

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
qGL_orig.total = zeros(numel(data),2); qOB_orig.total = zeros(numel(data),2);
qGL_calv.total = zeros(numel(data),2); qOB_calv.total = zeros(numel(data),2);
qGL_dhIS.total = zeros(numel(data),2); qOB_dhIS.total = zeros(numel(data),2);
qGL_dh.total = zeros(numel(data),2); qOB_dh.total = zeros(numel(data),2);
misfit = zeros(numel(data),2);
gsA = zeros(numel(data),1);
gsC = zeros(numel(data),1);
gaA = zeros(numel(data),1);
gaC = zeros(numel(data),1);

qGL_calvdh.total = zeros(numel(data),2); qOB_calvdh.total = zeros(numel(data),2);
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
    MUA = MUA_new;
    save("Delta_u.mat", "Delta_u","MUA","gsA","gsC","gaA","gaC","m","n");
end

CtrlVar=Ua2D_DefaultParameters;
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
                    yaxislabel = "\Delta q_{GL} [Gt/yr]";
        
                case 'Delta_qOB'
        
                    qorig = qOB_orig.total(:,cc);
                    ydata1 = qOB_calv.total(:,cc)-qorig;
                    ydata2 = qOB_dhIS.total(:,cc)-qorig;
                    ydata3 = qOB_dh.total(:,cc)-qorig;
                    ydata4 = qOB_calvdh.total(:,cc)-qorig;
                    yaxislabel = "\Delta q_{OB} [Gt/yr]";
        
                otherwise
        
                    continue;
        
            end
        
            if cc==1
                marker='o';                        
            elseif cc==2
                marker='s';               
            end

            s((cc-1)*4+1)=scatter(xdata,ydata1,markersize(:,1),marker,"filled",'MarkerEdgeColor','none');
            s((cc-1)*4+2)=scatter(xdata,ydata2,markersize(:,1),marker,"filled",'MarkerEdgeColor','none');
            s((cc-1)*4+3)=scatter(xdata,ydata3,markersize(:,1),marker,"filled",'MarkerEdgeColor','none');
            s((cc-1)*4+4)=scatter(xdata,ydata4,markersize(:,1),marker,"filled",'MarkerEdgeColor','none');

        end

        colororder([newcolors;newcolors]);

        xlim([floor(min(xdata)) ceil(max(xdata))]);
        
        for ii=0:4*cycles_to_plot(end)-1
            s(ii+1).AlphaData = alphavalue(:,floor(ii/4)+1);
            s(ii+1).MarkerFaceAlpha='flat';
        end
        
        for ii=1:4
            s(ii)=plot(0,0,'o','color',newcolors(ii,:),'MarkerSize',10,'MarkerFaceColor',newcolors(ii,:));
        end

        if numel(cycles_to_plot)==2   
            s(9)=plot(0,0,'ok','markersize',10);
            s(10)=plot(0,0,'sk','markersize',10);
            legend([s(1:4) s(9) s(10)],{'Calving','Ice Shelf thickness','Ice thickness','Calving + Ice thickness','no spinup','with spinup'},...
            'NumColumns',2,'Location','northwest');
        else
            legend(s(1:4),{'Calving','Ice Shelf thickness','Ice thickness','Calving + Ice thickness'},...
            'NumColumns',2,'Location','northwest');
        end

        
        
        xlabel(parameter_to_plot); ylabel(yaxislabel);
        ax=gca;
        
        if ismember(diagnostic_to_plot,["gsA","gsC","gaA","gaC"])
            ax.XScale='log';
        end
        ax.YScale='log';
        grid on;
        box on;


    case "Delta_u"
        
        fields_to_plot = fields(Delta_u);

        for cc=cycles_to_plot

            figure(cc*999);
            tlo=tiledlayout(2,numel(fields_to_plot),'TileSpacing','tight','TileIndexing', 'columnmajor');

            for ff=1:numel(fields_to_plot)

                deltau_av = mean(Delta_u.(fields_to_plot{ff}).map(:,:,cc),1);
                deltau_std = std(Delta_u.(fields_to_plot{ff}).map(:,:,cc),1);

                nexttile;
                PlotNodalBasedQuantities_JDR(gca,MUA.connectivity,MUA.coordinates,deltau_av(:),CtrlVar);
                colormap(othercolor('RdYlBu8'));
                title(fields_to_plot{ff});
                axis tight; axis off;
                caxis([-1000 1000]);

                for ff=1:numel(fields_to_plot)
                     cb=colorbar; cb.Layout.Tile='east'; %cb.Ylabel.String="average";
                end

                nexttile;
                PlotNodalBasedQuantities_JDR(gca,MUA.connectivity,MUA.coordinates,deltau_std(:),CtrlVar);
                colormap(othercolor('RdYlBu8'));
                axis tight; axis off;
                caxis([-100 100]);

                for ff=1:numel(fields_to_plot)
                     cb=colorbar; cb.Layout.Tile='east'; %cb.Ylabel.String="std";
                end
            end

        end
    
       

    otherwise 

        error("unknown case")
       
end


