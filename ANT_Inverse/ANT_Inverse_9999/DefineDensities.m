function  [UserVar,rho,rhow,g]=DefineDensities(UserVar,CtrlVar,MUA,time,s,b,h,S,B)
    
fprintf('Start loading densities...');

filename_densityfield = UserVar.DesnityInterpolant;
filename_densityfields = erase(filename_densityfield,"GriddedInterpolants_");
filename_densityfields = erase(filename_densityfield,".mat");
filename_densityfields = filename_densityfields + "_mesh_Nodes" + string(MUA.Nnodes) + "_Nele" + string(MUA.Nele) + ".mat";

if exist(filename_desnityfield,"file")
    load(filename_geometryfields,"rho")
else
    load(UserVar.GeometryInterpolants,'Frho');
    rho = Frho(MUA.coordinates(:,1),MUA.coordinates(:,2));
    rho(rho<100)=100;
    rho(rho>917)=917;
    save(filename_geometryfields,"rho","-append");
    clear Fhro;
end
	
rhow=1027; 
            
g=9.81/1000;

fprintf('done.\n');
   
end
