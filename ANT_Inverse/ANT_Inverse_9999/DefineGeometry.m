function [UserVar,s,b,S,B,alpha]=DefineGeometry(UserVar,CtrlVar,MUA,time,FieldsToBeDefined)

persistent FB Fs Fb;

if nargin<5
    FieldsToBeDefined='sbSB';
end

alpha=0 ;
fprintf('Loading s, b, S and B \n');

load(UserVar.GeometryInterpolants,'FB','Fs','Fb');

x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);

if contains(FieldsToBeDefined,'S')
    S = 0*x;
end

if contains(FieldsToBeDefined,'B')
    B = FB(x,y);
    B = inpaint_nans(B,4);
end

if contains(FieldsToBeDefined,'s')
    s = Fs(x,y);
    s = inpaint_nans(s,4);
end

if contains(FieldsToBeDefined,'b')
    b = Fb(x,y);
    b = inpaint_nans(b,4);
end


