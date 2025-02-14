function producefiles

addpath("../ANT_HelperFunctions/");

% Create_Geometry_from_reference_and_ds(datetime(2019,06,01),1);
% Klear;
% Velinterpolantfile = "GriddedInterpolants_1996-2003_MeaSUREs_ITSLIVE_Velocities";
% Geominterpolantfile = "GriddedInterpolants_Geometry_01-Jun-2019";
% Create_ExtrudedFields_GriddedInterpolants(Velinterpolantfile,Geominterpolantfile,[],1,'-geom-');
% Klear;

Create_MeaSUREs_ITSLIVE_20xx_VelocityGriddedInterpolants(2010,0);
Klear
Create_MeaSUREs_ITSLIVE_20xx_VelocityGriddedInterpolants(2011,0);
Klear
Create_MeaSUREs_ITSLIVE_20xx_VelocityGriddedInterpolants(2012,0);
Klear
Create_MeaSUREs_ITSLIVE_20xx_VelocityGriddedInterpolants(2013,0);
Klear
Create_MeaSUREs_ITSLIVE_20xx_VelocityGriddedInterpolants(2014,0);
Klear
Create_MeaSUREs_ITSLIVE_20xx_VelocityGriddedInterpolants(2015,0);
Klear
Create_MeaSUREs_ITSLIVE_20xx_VelocityGriddedInterpolants(2016,0);
Klear
Create_MeaSUREs_ITSLIVE_20xx_VelocityGriddedInterpolants(2017,0);
Klear
Create_MeaSUREs_ITSLIVE_20xx_VelocityGriddedInterpolants(2018,0);
Klear


return
% 
% Velinterpolantfile = "GriddedInterpolants_1996-2003_MeaSUREs_ITSLIVE_Velocities";
% Geominterpolantfile = "GriddedInterpolants_Geometry_01-Jun-2000";
% Create_ExtrudedFields_GriddedInterpolants(Velinterpolantfile,Geominterpolantfile,[],1,'-v-geom-');
% Klear;
% 
Velinterpolantfile = "GriddedInterpolants_2009-2010_MeaSUREs_ITSLIVE_Velocities";
Geominterpolantfile = "GriddedInterpolants_Geometry_01-Jun-2009";
Create_ExtrudedFields_GriddedInterpolants(Velinterpolantfile,Geominterpolantfile,[],1,'-v-');
Klear;
%
Velinterpolantfile = "GriddedInterpolants_2014-2015_MeaSUREs_ITSLIVE_Velocities";
Geominterpolantfile = "GriddedInterpolants_Geometry_01-Jun-2014";
Create_ExtrudedFields_GriddedInterpolants(Velinterpolantfile,Geominterpolantfile,[],1,'-v-');
Klear;
%
Velinterpolantfile = "GriddedInterpolants_2019-2020_MeaSUREs_ITSLIVE_Velocities";
Geominterpolantfile = "GriddedInterpolants_Geometry_01-Jun-2020";
Create_ExtrudedFields_GriddedInterpolants(Velinterpolantfile,Geominterpolantfile,[],1,'-v-');
Klear;
% 
% Velinterpolantfile = "GriddedInterpolants_2010-2011_MeaSUREs_ITSLIVE_Velocities";
% Geominterpolantfile = "GriddedInterpolants_Geometry_01-Jun-2010";
% Create_ExtrudedFields_GriddedInterpolants(Velinterpolantfile,Geominterpolantfile,1);
% Klear;
% 
% Velinterpolantfile = "GriddedInterpolants_2011-2012_MeaSUREs_ITSLIVE_Velocities";
% Geominterpolantfile = "GriddedInterpolants_Geometry_01-Jun-2011";
% Create_ExtrudedFields_GriddedInterpolants(Velinterpolantfile,Geominterpolantfile,1);
% Klear;
% 
% Velinterpolantfile = "GriddedInterpolants_2012-2013_MeaSUREs_ITSLIVE_Velocities";
% Geominterpolantfile = "GriddedInterpolants_Geometry_01-Jun-2012";
% Create_ExtrudedFields_GriddedInterpolants(Velinterpolantfile,Geominterpolantfile,1);
% Klear;
% 
% Velinterpolantfile = "GriddedInterpolants_2013-2014_MeaSUREs_ITSLIVE_Velocities";
% Geominterpolantfile = "GriddedInterpolants_Geometry_01-Jun-2013";
% Create_ExtrudedFields_GriddedInterpolants(Velinterpolantfile,Geominterpolantfile,1);
% Klear;
% 
% Velinterpolantfile = "GriddedInterpolants_2014-2015_MeaSUREs_ITSLIVE_Velocities";
% Geominterpolantfile = "GriddedInterpolants_Geometry_01-Jun-2014";
% Create_ExtrudedFields_GriddedInterpolants(Velinterpolantfile,Geominterpolantfile,1);
% Klear;
% 
%Velinterpolantfile = "GriddedInterpolants_2015-2016_MeaSUREs_ITSLIVE_Velocities";
%Geominterpolantfile = "GriddedInterpolants_Geometry_01-Jun-2015";
%Create_ExtrudedFields_GriddedInterpolants(Velinterpolantfile,Geominterpolantfile,[],1,'-v-');
%Klear;
% 
% Velinterpolantfile = "GriddedInterpolants_2016-2017_MeaSUREs_ITSLIVE_Velocities";
% Geominterpolantfile = "GriddedInterpolants_Geometry_01-Jun-2016";
% Create_ExtrudedFields_GriddedInterpolants(Velinterpolantfile,Geominterpolantfile,1);
% Klear;
% 
% Velinterpolantfile = "GriddedInterpolants_2017-2018_MeaSUREs_ITSLIVE_Velocities";
% Geominterpolantfile = "GriddedInterpolants_Geometry_01-Jun-2017";
% Create_ExtrudedFields_GriddedInterpolants(Velinterpolantfile,Geominterpolantfile,1);
% Klear;

% Velinterpolantfile = "GriddedInterpolants_2019-2020_MeaSUREs_ITSLIVE_Velocities";
% Geominterpolantfile = "GriddedInterpolants_Geometry_01-Jun-2019";
% Create_ExtrudedFields_GriddedInterpolants(Velinterpolantfile,Geominterpolantfile,1);
% Klear;
% 
% Velinterpolantfile = "GriddedInterpolants_2020-2021_MeaSUREs_ITSLIVE_Velocities";
% Geominterpolantfile = "GriddedInterpolants_Geometry_01-Jun-2020";
% Create_ExtrudedFields_GriddedInterpolants(Velinterpolantfile,Geominterpolantfile,1);
% Klear;
% 
% Velinterpolantfile = "GriddedInterpolants_2021-2022_MeaSUREs_ITSLIVE_Velocities";
% Geominterpolantfile = "GriddedInterpolants_Geometry_01-Jun-2021";
% Create_ExtrudedFields_GriddedInterpolants(Velinterpolantfile,Geominterpolantfile,1);
% Klear;


