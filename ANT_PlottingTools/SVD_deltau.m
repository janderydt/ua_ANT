function [U,S,V] = SVD_deltau(data,doplots)

if nargin==0
    load("Delta_u");
    data=Delta_u.Calv_dh.map(:,:,1);
    %data=Delta_u.dhIS.map(:,:,1);
    doplots=1;
elseif nargin==1
    doplots=1;
end

% check dimensions of input data 
% rows: nodes, columns: experiments
if numel(size(data))==2
    nx = size(data,1); ny = size(data,2);
    if nx < ny % number of experiments is highly likely to be smaller than number of nodes
        data = data'; % transpose data
    end
else
    error("check dimensions of input data");
end
data(isnan(data))=0;

[U,S,V] = svd(data,'econ');


if doplots

    CtrlVar=Ua2D_DefaultParameters;
    CtrlVar.PlotXYscale = 1e3;

    cumulative_variance = cumsum(diag(S).^2)/sum(diag(S).^2);
    
    % figure(111), tlo1=tiledlayout(2,2,'TileSpacing','tight'); title(tlo1,'A');
    % nexttile;
    % imagesc(A), axis off; colormap(slanCM('YlGnBu')); cb1=colorbar;
    % title('Original');
    
    % figure(222), tlo2=tiledlayout(2,2,'TileSpacing','tight'); title(tlo2,'t = 1 year');
    % nexttile; 
    % PlotNodalBasedQuantities_JDR(gca,MUA.connectivity,MUA.coordinates,A(:,1),CtrlVar), axis equal, axis off, colormap(slanCM('YlGnBu')); cb2=colorbar(gca);
    % title('Original');
    % 
    % figure(333), tlo3=tiledlayout(2,2,'TileSpacing','tight'); title(tlo3,'t = 850 years');
    % nexttile; 
    % PlotNodalBasedQuantities_JDR(gca,MUA.connectivity,MUA.coordinates,data(:,end),CtrlVar), axis equal, axis off, colormap(slanCM('YlGnBu')); cb3=colorbar(gca);
    % title('Original');
    
    figure(444); tlo4=tiledlayout(2,5,'TileSpacing','tight');
    for ii=1:10
        nexttile;
        PlotNodalBasedQuantities_JDR(gca,MUA.connectivity,MUA.coordinates,U(:,ii),CtrlVar);
        colormap(othercolor('RdYlBu8'));
        variance =100*[cumulative_variance(1); cumulative_variance(2:end)-cumulative_variance(1:end-1)];
        title(['mode',num2str(ii),' (',num2str(variance(ii),'%2.2f'),'%)']); caxis([-0.05 0.05]);
        axis tight; axis off;
    
        %cb=colormap; cb.visible='off';
    end
    cb4=colorbar(gca);
    plotind = 2;
    for r = [10 50 100]
        Xapprox = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
        %Xapprox(Xapprox<0)=0; 
        %Xapprox(Xapprox>1)=1;
        %figure(111); nexttile, imagesc(Xapprox), axis off, colormap(slanCM('YlGnBu'));
        %figure(222); nexttile, PlotNodalBasedQuantities_JDR(gca,MUA.connectivity,MUA.coordinates,Xapprox(:,1),CtrlVar), axis equal, axis off, colormap(slanCM('YlGnBu'));
        %title(['r=',num2str(r,'%d'),', ',num2str(100*energy(r),'%2.2f'),'% cumulative energy']);
        %figure(333); nexttile, PlotNodalBasedQuantities_JDR(gca,MUA.connectivity,MUA.coordinates,Xapprox(:,end),CtrlVar), axis equal, axis off, colormap(slanCM('YlGnBu'));
        %title(['r=',num2str(r,'%d'),', ',num2str(100*energy(r),'%2.2f'),'% cumulative energy']);
        %plotind = plotind + 1;
    end
    
    cb1.Layout.Tile = 'east';
    cb2.Layout.Tile = 'east';
    cb3.Layout.Tile = 'east';
    cb4.Layout.Tile = 'east';
    
    %% singular values
    figure, subplot(1,2,1)
    semilogy(diag(S),'k','linewidth',2); grid on;
    xlabel('r');
    ylabel('singlar value, \sigma_r');
    set(gca,'fontsize',14);
    subplot(1,2,2)
    plot(cumulative_variance,'k','linewidth',2); grid on;
    xlabel('r');
    ylabel('cumulative energy');
    set(gca,'fontsize',14);
end