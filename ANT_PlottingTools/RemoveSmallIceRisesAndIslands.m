function B = RemoveSmallIceRisesAndIslands(B)

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
    if A>250
        Bxtmp = [Bxtmp(:);xtmp(:);NaN];
        Bytmp = [Bytmp(:);ytmp(:);NaN];
    end
end
B.x{1} = Bxtmp;
B.y{1} = Bytmp;
