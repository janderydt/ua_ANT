function [UserVar,s,b,S,B,alpha]=DefineGeometry(UserVar,CtrlVar,MUA,time,FieldsToBeDefined)

persistent FB Fs Fb;

if nargin<5
    FieldsToBeDefined='sbSB';
end

s=[]; b=[]; S=[]; B=[];
alpha=0 ;

fprintf('Loading geometry fields \n');

x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);

if contains(FieldsToBeDefined,'S')
    S = 0*x;
end

if contains(FieldsToBeDefined,'B')
    if isempty(FB)
        load(UserVar.GeometryInterpolants,'FB');
    end
    B = FB(x,y);
    B = inpaint_nans(B,4);
end

if contains(FieldsToBeDefined,'s')
    if isempty(Fs)
        load(UserVar.GeometryInterpolants,'Fs');
    end
    s = Fs(x,y);
    s = inpaint_nans(s,4);
end

if contains(FieldsToBeDefined,'b')
    if isempty(Fb)
        load(UserVar.GeometryInterpolants,'Fb');
    end
    b = Fb(x,y);
    b = inpaint_nans(b,4);
end


