function UserVar=UaOutputs(UserVar,CtrlVar,MUA,BCs,F,l,GF,InvStartValues,InvFinalValues,Priors,Meas,BCsAdjoint,RunInfo)

v2struct(F);
time=CtrlVar.time; 
    
% check if folder 'UaOutputsDirectory' exists, if not create
if exist(fullfile(cd,UserVar.UaOutputDirectory),'dir')~=7
    mkdir(UserVar.UaOutputDirectory) ;
end

% matlab output
YearsSinceStart = num2str(sprintf('%04d',time));
FileName=sprintf([UserVar.UaOutputDirectory,'/',CtrlVar.Experiment,'_',YearsSinceStart]);
fprintf(' Saving data in %s \n',FileName);
save(FileName,'UserVar','CtrlVar','MUA','time','s','b','S','B','h','ub','vb','C','dhdt','AGlen','m','n','rho','rhow','as','ab','GF');

% netcdf output
% ncfilename = [UserVar.UaOutputDirectory,'/',CtrlVar.Experiment,'_',YearsSinceStart,'.nc'];
% if exist(ncfilename,'file')
%     delete(ncfilename);
% end
% dimensions.name={'nx'};
% dimensions.dim=[MUA.Nnodes];
% [lat,lon] = psxy2ll(MUA.coordinates(:,1),MUA.coordinates(:,2),-71,0);
% variables.name={'lon','lat','iceThickness','upperIceSurface','lowerIceSurface','bedElevation','seaSurface'};
% variables.dim=[1; 1; 1; 1; 1; 1; 1];
% variables.units={'m','m','m','m','m','m','m'};
% variables.description={'x location of Ua nodes','y location of Ua nodes','ice thickness',...
%     'upper surface elevation above geoid','lower surface elevation above geoid','bed elevation above geoid','sea surface elevation'};
% 
% for ii=1:length(variables.name)
%     I = find(variables.dim(ii,:)); 
%     dim = {dimensions.name{I(1)},dimensions.dim(I(1))};
%     for jj=2:length(I)
%         dim = {dim{:}, dimensions.name{I(jj)},dimensions.dim(I(jj))};
%     end
%     nccreate(ncfilename,variables.name{ii},...
%          'Dimensions',dim,'FillValue', NaN,...
%          'Format','netcdf4','Datatype','single');
%     ncwriteatt(ncfilename,variables.name{ii},'units',variables.units{ii});
%     ncwriteatt(ncfilename,variables.name{ii},'description',variables.description{ii});
% end
% ncwriteatt(ncfilename,'/','_FillValue',NaN);
% ncwrite(ncfilename,'lon',lon);
% ncwrite(ncfilename,'lat',lat);
% ncwrite(ncfilename,'iceThickness',h);
% ncwrite(ncfilename,'upperIceSurface',s);
% ncwrite(ncfilename,'lowerIceSurface',b);
% ncwrite(ncfilename,'bedElevation',B);
% ncwrite(ncfilename,'seaSurface',S);

