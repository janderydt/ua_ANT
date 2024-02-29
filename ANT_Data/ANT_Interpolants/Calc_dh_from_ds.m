function [b,h] = Calc_dh_from_ds(X,Y,s,B,rho,rhow)

S = 0*X;

% get a rough and a reasonable initial estimate for b if none is provided
% The lower surface b is 
%
%
%   b=max( B , (rhow S - rho s)/(rhow-rho) ) 
%   where
%
%  h_f = rhow (S-B) / rho
%
%  b=s-h_f =

hf=rhow*(S-B)./rho ;

b0 =  max(B,(rho.*s-rhow.*S)./(rho-rhow)) ; 

b=b0;
h=s-b;

% iteration
ItMax=30 ; tol=100*eps ;  J=Inf ; I=0 ;
JVector=zeros(ItMax,1)+NaN ;        

CtrlVar.kH = 1;
CtrlVar.Hh0 = 0;

while I< ItMax && J > tol
    I=I+1;

    G = HeavisideApprox(CtrlVar.kH,h-hf,CtrlVar.Hh0);  % 1
    dGdb=-DiracDelta(CtrlVar.kH,h-hf,CtrlVar.Hh0) ;
    
    F0=    b - G.*B - (1-G).*(rho.*s-rhow.*S)./(rho-rhow) ;
    dFdb = 1 - dGdb.* (B -  (rho.*s-rhow.*S)./(rho-rhow)) ;
    
    db= -F0./dFdb ;
    
    b=b+db ;
    h=s-b ;
    
    F1 =    b - G.*B - (1-G).*(rho.*s-rhow.*S)./(rho-rhow) ;
    
    JLast=J ;
    J=sum(F1.^2,'all','omitnan')/2 ;
    
    JVector(I)=J ;
    
end

fprintf("%s iterations needed for bh calculation.\n",string(I));

end


