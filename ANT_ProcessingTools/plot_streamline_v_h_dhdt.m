function plot_streamline_v_h_dhdt

persistent SL Fv Fh FGF

addpath(getenv("froot_ua")+"cases/ANT");
addpath(getenv("froot_tools"));

glaciers = "PIG";%["PIG","Thwaites","Pope","SmithEast","SmithWest","Kohler"];

UserVar.home = "/mnt/md0/Ua/cases/ANT/";
UserVar.type = "Diagnostic";

domain = "AMUND";
years = ["2000-2020"]; % can be individual years, e.g. "2000", or multiple years
    % seperated by a "-", e.g. "2000-2020"
cycle = [1];  % specify a cycle for each string of years above
slidinglaw = ["Weertman"]; % specify a sliding law for each string of years above

xmin = zeros(numel(glaciers),1)+200;
xmax = zeros(numel(glaciers),1);

for cc=1:numel(years)

    load("inversiondata_"+domain+"_"+slidinglaw(cc)+".mat");
    data_inv = data;
    
    individual_years = strsplit(years(cc),"-");

    for tt=1:numel(individual_years)

        year_tmp = individual_years(tt);
        if year_tmp == "2000"
            UserVar.Table = UserVar.home+"ANT_Diagnostic/RunTable_ARCHER2_Diagnostic_"+domain+"_Original_CalvThick_"+...
                slidinglaw(cc)+"_2020.csv";
        else
            UserVar.Table = UserVar.home+"ANT_Diagnostic/RunTable_ARCHER2_Diagnostic_"+domain+"_Original_CalvThick_"+...
                slidinglaw(cc)+"_"+year_tmp+".csv";
        end
    
        % read run table
        RunTable = ANT_ReadWritetable(UserVar,UserVar.Table,[],'read');
        
        % ExpIDs
        ExpID = RunTable{:,"ExpID"};
        idrange = [min(ExpID) max(ExpID)];
        Ind = find(ExpID>=idrange(1) & ExpID<=idrange(2));
        
        % only keep experiments that have finished
        Ind_finished = RunTable{Ind,"Finished"}==1;
        Comments = RunTable{Ind,"Comments"};
        Cycle = RunTable{Ind,"InverseCycleA"};
        
        if year_tmp == "2000"
            Ind_Orig = contains(Comments,"Original geometry 2000");
            Ind = Ind(Ind_finished==1 & Ind_Orig==1 & Cycle == cycle(cc));
        else
            Ind_CalvThick = contains(Comments,"Ice front and thickness "+year_tmp);
            Ind = Ind(Ind_finished==1 & Ind_CalvThick==1 & Cycle == cycle(cc));
        end
    
        kk=1;
        
        %% Gather data
        for ii=1:numel(Ind)
    
            inverse_experiment = RunTable{Ind(ii),"InverseA"};
            Ind_inv = find([data_inv(:).InverseExpID]==inverse_experiment);
            m(kk) = data_inv(Ind_inv).m;
            n(kk) = data_inv(Ind_inv).n;
            gaA(kk) = data_inv(Ind_inv).gaA;
            gsA(kk) = data_inv(Ind_inv).gsA;
            gaC(kk) = data_inv(Ind_inv).gaC;
            gsC(kk) = data_inv(Ind_inv).gsC;
    
            folder = UserVar.home+"/ANT_Diagnostic/cases/"+domain+"_nsmbl_Diagnostic_"+ExpID(Ind(ii));
            outputfiles = dir(folder+"/ResultsFiles/*.mat");
    
            if ~isempty(outputfiles)
    
                outputfile = outputfiles(1).folder + "/" + outputfiles(1).name;
                expinfo = Comments{ii};
    
                if isempty(SL)
                    n=5;
                    SL = Define_Streamline(outputfile,glaciers);
                    load(outputfile,"MUA","F");
                    x = MUA.coordinates(:,1); y = MUA.coordinates(:,2);
                    Fv = scatteredInterpolant(x,y,hypot(F.ub,F.vb));
                    Fh = Fv; Fh.Values = F.h;
                    FGF = Fv; FGF.Values = F.GF.node;
                else
                    load(outputfile,"F");
    
                    Fv.Values = hypot(F.ub,F.vb);
                    Fh.Values = F.h;
                    FGF.Values = F.GF.node;
                end
                
                %MUA=UpdateMUA(CtrlVarInRestartFile,MUA);
                %[~,dhdt]=dhdtExplicitSUPG(UserVarInRestartFile,CtrlVarInRestartFile,MUA,F,BCs);
                %Fdhdt = scatteredInterpolant(x,y,dhdt);
    
                for gg=1:numel(glaciers)
                    SL(gg).V(kk,:) = Fv(SL(gg).x,SL(gg).y);
                    SL(gg).H(kk,:) = Fh(SL(gg).x,SL(gg).y);
                    %dHdT = Fdhdt(SLx,SLy);
                    GF = FGF(SL(gg).x,SL(gg).y);
                    I = find(GF<0.9); 
                    SL(gg).SLd_GL(kk) = SL(gg).d(I(1));
                end
    
                kk=kk+1;
    
            end
    
            disp("Done "+string(ii)+" out of "+string(numel(Ind)));
    
        end  
    
        %% load observations
        addpath(getenv("froot_data")+"Measures/Measures_annual");
        
        if year_tmp=="2000"
            fname = "GriddedInterpolants_1996-2003_MeaSUREs_ITSLIVE_Velocities.mat";
        else
            fname = "GriddedInterpolants_"+string(double(year_tmp)-1)+"-"+string(double(year_tmp))+...
                "_MeaSUREs_ITSLIVE_Velocities.mat";
        end
        load("../ANT_Data/ANT_Interpolants/"+fname);
        v_tmp = hypot(Fus.Values,Fvs.Values);
        Fu = Fus; Fu.Values = v_tmp;
        std_tmp = hypot(Fxerr.Values,Fyerr.Values);
        Fstd = Fus; Fstd.Values = std_tmp;
        for gg=1:numel(glaciers)
            SL(gg).vMeas = Fu(SL(gg).x,SL(gg).y);
            SL(gg).stdMeas = Fstd(SL(gg).x,SL(gg).y);
        end
        
        %% set up figures
        nrows = numel(years);
        if any(cycle>1)
            nrows = nrows+1;
        end
        ncolumns = numel(individual_years);

        for gg=1:numel(glaciers)
            
            if cc==1 && tt==1
                Hfig(100+gg)=fig('units','inches','width',120*12/72.27,'height',65*12/72.27,'fontsize',14,'font','Helvetica');
           
                tlo_fig(100+gg) = tiledlayout(Hfig(100+gg),nrows,ncolumns,"TileSpacing","compact");
                for i = 1:(nrows)*ncolumns
                    ax_fig(100+gg,i) = nexttile(tlo_fig(100+gg),i); hold on;
                end
            end
            
            for ii=1:size(SL(gg).V,1)
    
                indmin = min(m); indmax = max(m);
                colorind = round(1+(64-1)/(indmax-indmin)*(m(ii)-indmin));
    
                plot(ax_fig(100+gg,ncolumns*(cc-1)+tt),SL(gg).d/1e3,SL(gg).V(ii,2:end),LineStyle='-',LineWidth=1,Color=[1 0.5 0.5]);%Color=CM(colorind,:));
                plot(ax_fig(100+gg,ncolumns*(cc-1)+tt),[SL(gg).SLd_GL(ii),SL(gg).SLd_GL(ii)]/1e3,[0 10000],'-k');
        
                %plot(ax_fig(2),SL(gg).d/1e3,SL(gg).H(ii,2:end),LineStyle=style,LineWidth=width,Color=CM(ceil(ii/2),:));
                %plot(ax_fig(2),[SL(gg).SLd_GL(ii),SL(gg).SLd_GL(ii)]/1e3,[0 2000],'-k');
            end
        
            plot(ax_fig(100+gg,ncolumns*(cc-1)+tt),SL(gg).d/1e3,SL(gg).vMeas(2:end),LineStyle="-",LineWidth=2,Color="k");
            plot(ax_fig(100+gg,ncolumns*(cc-1)+tt),SL(gg).d/1e3,SL(gg).vMeas(2:end)-SL(gg).stdMeas(2:end),LineStyle="--",LineWidth=2,Color="k");
            plot(ax_fig(100+gg,ncolumns*(cc-1)+tt),SL(gg).d/1e3,SL(gg).vMeas(2:end)+SL(gg).stdMeas(2:end),LineStyle="--",LineWidth=2,Color="k");
        
            xlim(ax_fig(100+gg,ncolumns*(cc-1)+tt),[0 150]);
            if tt==1
                ylim(ax_fig(100+gg,ncolumns*(cc-1)+tt),[0 max(SL(gg).V(:))]);
                ylabel(ax_fig(100+gg,ncolumns*(cc-1)+tt),"Flowline speed [m/yr]");
            else
                ylims = ylim(ax_fig(100+gg,1));
                ylim(ax_fig(100+gg,ncolumns*(cc-1)+tt),ylims);
            end
            grid(ax_fig(100+gg,ncolumns*(cc-1)+tt),"on");
            box(ax_fig(100+gg,ncolumns*(cc-1)+tt),"on");
            title(ax_fig(100+gg,ncolumns*(cc-1)+tt),year_tmp+" (cycle "+string(cycle(cc))+", "+...
                slidinglaw(cc)+")");

            if any(cycle>1)
                edges = linspace(min(SL(gg).d/1e3),max(SL(gg).d/1e3),40);
                xmid = 0.5*(edges(1:end-1)+edges(2:end));
                wbar = xmid(2)-xmid(1);
                [n,~,~] = histcounts(SL(gg).SLd_GL(:)/1e3,edges,'Normalization','percentage'); 
                bar(ax_fig(100+gg,(nrows-1)*ncolumns+tt),xmid+(-1)^cc*wbar/4,n,0.4,'FaceAlpha',0.5);
                xlim(ax_fig(100+gg,(nrows-1)*ncolumns+tt),[0 150]);
                ylim(ax_fig(100+gg,(nrows-1)*ncolumns+tt),[0 100]);
                grid(ax_fig(100+gg,(nrows-1)*ncolumns+tt),"on");
                box(ax_fig(100+gg,(nrows-1)*ncolumns+tt),"on");
                ylabel(ax_fig(100+gg,(nrows-1)*ncolumns+tt),"GL location percentage");
            end
    
            xlabel(tlo_fig(100+gg),'Distance [km]');
            title(tlo_fig(100+gg),glaciers(gg));

        end

        SL=[];
    
    end

end

end

%% set up streamline data
function SL = Define_Streamline(file,glaciers)

load(file);
x=MUA.coordinates(:,1); y=MUA.coordinates(:,2); xmin=min(x); xmax=max(x); ymin=min(y); ymax=max(y);
[X,Y]=meshgrid([xmin:1e3:xmax],[ymin:1e3:ymax]);
Fu=scatteredInterpolant(x,y,F.ub);
Fv = Fu; Fv.Values = F.vb;
U=Fu(X,Y);
V=Fv(X,Y);
% remove nan from Boundary
I = find(isnan(MUA.Boundary.x));
if ~isempty(I)
    Bx = MUA.Boundary.x(1:I(1)-1);
    By = MUA.Boundary.y(1:I(1)-1);
else
    Bx = MUA.Boundary.x;
    By = MUA.Boundary.y;
end
U(find(~inpoly2([X(:),Y(:)],[Bx(:) By(:)])))=nan;
V(find(~inpoly2([X(:),Y(:)],[Bx(:) By(:)])))=nan;

FigSL=figure; hold on;

for ii=1:numel(glaciers)
    switch glaciers(ii)
        case "PIG"
            xstart = -1581780;
            ystart = -168332;
        case "Thwaites"
            xstart = -1457650;
            ystart = -444298;
        case "Pope"
            xstart = -1449590;
            ystart = -574470;
        case "SmithEast"
            xstart = -1465410;
            ystart = -617659;
        case "SmithWest"
            xstart = -1476460;
            ystart = -643130;
        case "Kohler"
            xstart = -1468270;
            ystart = -664000;
    end
 
    SL_tmp=streamline(gca(FigSL),X,Y,U,V,xstart,ystart); hold on;
    axis equal; axis tight;
    plot(Bx,By,'-k');
    PlotGroundingLines(CtrlVar,MUA,F.GF,[],[],[],'-k','linewidth',1);

    SLx=SL_tmp.XData; SLy=SL_tmp.YData;
    % remove SL outside Ua boundary
    I = find(~inpoly2([SLx(:) SLy(:)],[Bx(:) By(:)]));
    SLx(I)=[]; SLy(I)=[];
    SLd=cumsum(hypot(SLx(2:end)-SLx(1:end-1),SLy(2:end)-SLy(1:end-1)));

    SL(ii).x = SLx;
    SL(ii).y = SLy;
    SL(ii).d = SLd;
end

end



