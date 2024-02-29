function B = Define_UaNodes_PerBasin(MUA,B)

%% detect basin for each node
% 1. inpoly
MUA_tmp = MUA;
for ii=1:numel(B.x)
    Bx = B.x{ii};
    By = B.y{ii};
    IB(ii).ele = [];
    Inan = [0; find(isnan(Bx)); numel(Bx)];
    for jj=1:numel(Inan)-2
        IB(ii).ele = [IB(ii).ele; find(inpoly([MUA_tmp.xEle(:) MUA_tmp.yEle(:)],...
            [Bx(Inan(jj)+1:Inan(jj+1)-1) By(Inan(jj)+1:Inan(jj+1)-1)]))];
    end
    % avoid counting double
    MUA_tmp.xEle(IB(ii).ele,:) = 1e12;
    MUA_tmp.yEle(IB(ii).ele,:) = 1e12;
end
% 2. any remaining nodes are allocated to basin with nearest boundary
remaining = find(MUA_tmp.xEle~=1e12);
PQ = [MUA_tmp.xEle(remaining) MUA_tmp.yEle(remaining)]; dist = [];
for ii=1:numel(B.x)
    P = [B.x{ii} B.y{ii}];
    [~,dist(ii,:)] = dsearchn(P,PQ);
end
[~,I] = min(dist,[],1);
    
for ii=1:numel(B.x)
    IB(ii).ele = [IB(ii).ele; remaining(find(I==ii))];
    B.UaEle{ii} = IB(ii).ele;   
end
