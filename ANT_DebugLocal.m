function ANT_DebugLocal

% script that allows you to continue an existing simulation in the GUI
% on your local machine. This is mostly useful for debugging purposes.

UserVar.home = pwd;
UserVar.type = "Inverse";
UserVar.ExpID = 10774;
UserVar.domain = "ANT_nsmbl";

RunFolder = UserVar.home+"/ANT_"+UserVar.type+"/cases/"+UserVar.domain+"_"+UserVar.type+...
    "_"+string(UserVar.ExpID);
files = dir(RunFolder+"/RunTable*.csv");
TableToRead = files(1).folder+"/"+files(1).name;

% look for ExpID in runtable
RunTable = ANT_ReadWritetable(UserVar,TableToRead,[],'read');

% adjust runtable
RunTable{1,"Submitted"}=0;
RunTable{1,"Running"}=0;
RunTable{1,"Error"}=0;
RunTable{1,"Restart"}=1;

TableToWrite = UserVar.home+"/ANT_"+UserVar.type+"/RunTable_debug.csv";
ANT_ReadWritetable(UserVar,TableToWrite,RunTable,'write');

addpath(getenv("froot_ua")+"/cases/ANT/");
addpath(getenv("froot_ua")+"/cases/ANT/ANT_"+UserVar.type);

cd(UserVar.home+"/ANT_"+UserVar.type);
ANT_UaWrapper("ua_config_debug.txt",9999,UserVar.type,1,UserVar.ExpID);