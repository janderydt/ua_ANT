function B = RemoveIceRisesAndIslands(B,cutoff)

% cutoff (in km2) is used to remove ice rises and islands with ab area
% larger than the cutoff

if nargin < 2
    cutoff = 250; %remove small ice rises and islands.
end

% Only keep ice rises and islands larger than 250km2. We don't treat
% smaller ice rises and/or islands seperately
Bx = B.x{1}; Bxtmp = [];
By = B.y{1}; Bytmp = [];
I = [0; find(isnan(Bx)); numel(Bx)+1];
for ii=1:numel(I)-1
    xtmp = Bx(I(ii)+1:I(ii+1)-1);
    ytmp = By(I(ii)+1:I(ii+1)-1);
    P = polyshape(xtmp(:),ytmp(:));
    A = area(P)/1e6; % square km
    if A>cutoff
        Bxtmp = [Bxtmp(:);xtmp(:);NaN];
        Bytmp = [Bytmp(:);ytmp(:);NaN];
    end
end
if isempty(Bxtmp)
    B.area_km2(1)=[];
    B.name=B.name(2:19);
    B.x_center(1)=[];
    B.x=B.x(2:19);
    B.y_center(1)=[];
    B.y=B.y(2:19);
else
    B.x{1} = Bxtmp;
    B.y{1} = Bytmp;
end
