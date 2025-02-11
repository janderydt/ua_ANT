function  BCs=DefineBoundaryConditions(UserVar,CtrlVar,MUA,BCs,time,s,b,h,S,B,ub,vb,ud,vd,GF)

  
if ~isempty(UserVar.BaseMesh.FixedBoundaryPoints)
         
    % load points that define the line-segments along which the BCs are to be defined
    load(UserVar.BaseMesh.FixedBoundaryPoints);
    AA=[xOuter(1:end-1) yOuter(1:end-1)] ; BB=[xOuter(2:end) yOuter(2:end)];

    % find all boundary nodes within 1m distance from the line segment.
    tolerance=1;
    I = DistanceToLineSegment([x(MUA.Boundary.Nodes) y(MUA.Boundary.Nodes)],AA,BB,tolerance);
    
    BCs.vbFixedNode=MUA.Boundary.Nodes(I);
    BCs.ubFixedNode=MUA.Boundary.Nodes(I);
    
    load(UserVar.VelocityInterpolants);

    BCs.ubFixedValue=Fus(x(BCs.ubFixedNode),y(BCs.ubFixedNode));
    BCs.vbFixedValue=Fvs(x(BCs.vbFixedNode),y(BCs.vbFixedNode));
    nanInd = find(isnan(BCs.ubFixedValue) | isnan(BCs.vbFixedValue));
    BCs.ubFixedValue(nanInd)=0;
    BCs.vbFixedValue(nanInd)=0;

end

if UserVar.SpinupCycle && UserVar.Spinup.Cycle > 1

    %% spinup: keep ice thickness for ice shelves fixed
    if CtrlVar.AdaptMesh==0

        IndIS= find(GF.node < 0.5);
        BCs.hFixedNode=IndIS;
        BCs.hFixedValue=h(IndIS);

        fprintf("Adding fixed ice shelf thickness condition to BCs.\n");

    else
    
        fprintf("Implement fixed ice shelf thickness for a changing mesh.\n");
        error("Implement fixed ice shelf thickness for a changing mesh");
        
    end


end

end