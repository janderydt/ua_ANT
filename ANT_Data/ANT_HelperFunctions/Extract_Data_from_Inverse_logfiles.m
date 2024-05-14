function Extract_Data_from_Inverse_logfiles

%% setup Ua folders
[~,hostname]= system("hostname"); 
if strfind(hostname,"C23000100")
    run /mnt/md0/Ua/setup_Ua2D.m
    UserVar.hostname = "C23000100";
    NumWorkers = 25;
end
UserVar.Table = "/home/wchm8/Documents/Ua/cases/ANT/ANT_Inverse/RunTable.csv";
UserVar.type = "Inverse";

addpath(genpath(getenv("froot_ua")+"/cases/ANT/"));

%% read table
RunTable = ANT_ReadWritetable(UserVar,[],'read');

Iexisting = find(RunTable{:,'ExpID'}~=0);

for ii=1:numel(Iexisting)

    ind = Iexisting(ii);

    ExpName = [char(RunTable{ind,'Domain'}),'_',char(UserVar.type),'_',num2str(RunTable{ind,'ExpID'})];

    x = fileread(getenv("froot_ua")+"/cases/ANT/ANT_"+UserVar.type+"/"+ExpName+"/"+ExpName+".out");
    xc = strsplit(x,sprintf('\n')); % split to cellstr
    xc_strip = strip(xc); % remove leading and trailing white spaces
    ix = regexp(xc_strip,'^\d'); % match a digit at beginning of line
    idxdigit = ~cellfun(@isempty,ix); % which are digits?
    iterations=xc_strip(idxdigit); % extract";
    last_iteration = iterations{end};

    fprintf("%s: %s.\n",ExpName,string(last_iteration));
       
end